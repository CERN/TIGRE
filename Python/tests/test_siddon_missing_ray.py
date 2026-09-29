# -*- coding: utf-8 -*-
"""Siddon forward projection must TERMINATE when rays miss the volume.

Regression test for a hang in `Common/CUDA/Siddon_projection.cu`, which is
compiled into both the Python extension and the MATLAB mex (see
`MATLAB/Compile.m`), so this covers both bindings.

THE BUG

    // line intersects voxel space ->   am<aM
    if (am>=aM)
        detector[idx]=0;        // <-- no return; execution falls through

`am>=aM` means the ray never enters the volume. Without a `return`,
execution continues to

    unsigned long Np=(imax-imin+1)+(jmax-jmin+1)+(kmax-kmin+1);
    for (unsigned long ii=0; ii<Np; ii++){ ... }

where imin..kmax are only meaningful for a ray that does enter. Measured on
the ray that hangs (proj 53, row 619, col 313 of the geometry below) by
instrumenting the kernel to dump its intermediates:

  1. the eps guard snaps ray.x = 0.00025 and ray.z = -4e-6 to exactly 0;
     the ray lies along -Y to ~1e-6 rad;
  2. source.x = 280.4 is OUTSIDE nVoxelX = 64, so -source.x/0 and
     (nVoxelX-source.x)/0 are BOTH -inf; fminf and fmaxf collapse onto one
     value and axm == axM == -inf;
  3. aM = fminf(axM, ayM, azM) = -inf while am = 0.51 -- the miss IS
     detected correctly;
  4. without the return the Z block runs anyway:
     kmin = ceilf(source.z + aM*ray.z), and (-inf)*0 = NaN;
  5. Np = NaN, and on CUDA (unsigned long)NaN = 9223372036854775808 = 2**63.
     That thread iterates ~9e18 times; the kernel never returns.

There is no error, no warning and no result. On the GPU it shows as 100% SM
occupancy with an IDLE memory controller and no PCIe traffic, because one or
two threads are re-fetching the same clamped texture address while every
other thread finished long ago. It is 2 rays in 261,479,799 here.

A NEGATIVE Np does NOT do this. Measured: CUDA saturates negatives and -inf
to 0, a zero-trip loop. An earlier version of this text said the sum wrapped
to ~1.8e19; that was wrong. It is specifically NaN.

Figure: Frontispiece/siddon_missed_ray_hang.png (generator alongside) draws
this ray in the kernel's voxel coordinates with the real captured numbers.

WHEN IT TRIPS

One ray needs all of these at once:

  * TWO direction components under eps -- the ray lies along a voxel axis.
    One supplies the +-inf, the other multiplies it: inf*0 = NaN.
  * the source OUTSIDE the volume's extent along a zeroed axis. Inside, the
    two quotients are -inf and +inf, stay apart, and are harmless.
  * the ray MISSES the volume, so the fall-through path is taken at all.

Why both `nDetector` entries odd: an odd count puts a pixel EXACTLY on the
detector's centre line, and when that aligns with the source's z the whole
centre row has ray.z == 0 -- 1179 columns x 179 views = 211,000 rays here.
Both hanging rays sat on row 619 = (1239-1)/2. An even count straddles the
centre by half a pixel (0.27 voxels, 270x eps) and never snaps. The second
zero is then an angular coincidence somewhere along that row; adding one
pixel to either dimension removes it, removing one re-creates it.

Why the object offset matters: the geometry below is a 100x125x40 mm volume
with DSO ~400-500 mm. Centred (offset 0) it does not hang; offset
+-50/100/200 mm it hangs; +-300/433 mm it does not hang again. It is a BAND
-- where the coincidence happens to land for this sampling -- not a
threshold.

Circular CT is not immune in principle: an axis-parallel ray from a source
at radius DSO also misses. But off-axis slab geometry with a large detector
sweeps ray directions widely enough that the ~1e-6 rad coincidence actually
lands, which is why it surfaces in tomosynthesis and laminography, and why a
stock centred cone-beam geometry with an odd x odd detector did NOT
reproduce it here.

THE FIX

    if (am>=aM){
        detector[idx]=0;
        return;
    }

Verified value-neutral: every configuration that returned before the fix
returns bit-identical sums after it.

THE PARALLEL KERNEL IS NOT COVERED HERE, AND NOT FOR WANT OF TRYING

`Siddon_projection_parallel.cu` carries the identical hazard --

    unsigned long Np=(imax-imin+1)+(jmax-jmin+1);

-- and was given the same `return`, but the hang could NOT be reproduced in
it. 26 configurations on a pre-fix build all returned in ~0.75 s: both
detector parities, offOrigin from +-5 mm to +-2000 mm, spanning the band
that hangs the cone kernel and well past it. Dispatch was confirmed genuine
(`_Ax.pyx` sends mode=="parallel" to siddon_ray_projection_parallel, and the
sums differ from cone mode), so the parallel path really was exercised.

Two plausible reasons it resists, neither established:
  * only TWO terms sum into Np instead of three, so there is less room to
    go negative;
  * that kernel already carries defensive code the cone one lacks -- an
    epsilon on the divisions and explicit copysignf(1e11,...) handling for
    ray.x/ray.y == 0 -- added, per its own comment, because "in paralel
    beam, often ray.y or ray.x=0; this leads to infinities propagating and
    breaking everything". It has been hardened once already.

So the parallel fix is justified BY INSPECTION and is UNVERIFIED. Do not
read the cone tests below as covering it. If you find a parallel trigger,
add it here; the search above is recorded so it need not be repeated.

Each projection runs in a CHILD PROCESS with a timeout, so a regression
FAILS this test rather than hanging the suite.
"""
import multiprocessing as mp

import numpy as np

TIMEOUT_S = 90

# MUST match the configuration the hang was verified on. The offsets and the
# per-view DSD/DSO are interpolated across this many views, so changing it
# resamples the trajectory and the grazing view that trips the bug may not
# be visited at all: at 60 views this test PASSES on a known-broken build.
# Do not lower it to make the test faster.
N_ANGLES = 179


def _geometry(n_det, off_mm):
    import tigre
    geo = tigre.geometry(mode="cone", default=True,
                         nVoxel=np.array([128, 128, 64]))
    geo.sVoxel = np.array([100.0, 125.0, 40.0])
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    geo.nDetector = np.asarray(n_det)
    geo.dDetector = np.array([0.42, 0.42])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.DSD = np.linspace(580.0, 920.0, N_ANGLES)
    geo.DSO = np.linspace(407.0, 514.0, N_ANGLES)
    geo.offOrigin = np.zeros((N_ANGLES, 3))
    geo.offOrigin[:, 2] = np.linspace(-off_mm, off_mm, N_ANGLES)
    geo.offDetector = np.zeros((N_ANGLES, 2))
    geo.offDetector[:, 1] = np.linspace(-1.7 * off_mm, 1.7 * off_mm, N_ANGLES)
    angles = np.linspace(0.9163, 2.2253, N_ANGLES).astype(np.float32)
    return geo, angles


def _project(n_det, off_mm, q):
    import tigre
    geo, angles = _geometry(n_det, off_mm)
    img = np.zeros(tuple(int(v) for v in geo.nVoxel), dtype=np.float32)
    img[32:96, 32:96, 16:48] = 1.0
    proj = tigre.Ax(img, geo, angles)
    q.put((float(np.asarray(proj).sum()),
           bool(np.all(np.isfinite(np.asarray(proj))))))


def _run(n_det, off_mm):
    """Project in a child process. Returns (sum, all_finite) or None on
    timeout -- None IS the failure this test exists to catch."""
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=_project, args=(n_det, off_mm, q))
    p.start()
    p.join(TIMEOUT_S)
    if p.is_alive():
        p.terminate()
        p.join()
        return None
    return None if q.empty() else q.get()


def test_missed_rays_terminate_with_odd_detector():
    """The reported case: both detector dimensions odd, object off-axis."""
    for off_mm in (50.0, 100.0, 200.0):
        got = _run(np.array([1239, 1179]), off_mm)
        assert got is not None, (
            f"Ax did not return within {TIMEOUT_S}s at offOrigin +-{off_mm} mm "
            f"with a 1239x1179 detector -- the missed-ray fall-through in "
            f"Siddon_projection.cu has regressed")
        total, finite = got
        assert finite, f"non-finite projection values at +-{off_mm} mm"
        assert total > 0.0, f"empty projection at +-{off_mm} mm: {total}"


def test_detector_parity_does_not_change_the_result():
    """Padding the detector by one pixel must not change what is measured.

    This is the other half: the fix has to make the hang go away WITHOUT
    moving any value. Odd and even detectors are different samplings, so
    the totals differ slightly -- but only by the added row/column, not by
    a factor.
    """
    ref = _run(np.array([1240, 1180]), 100.0)
    assert ref is not None, "even x even detector did not return"
    for n_det in ([1239, 1179], [1239, 1180], [1240, 1179]):
        got = _run(np.array(n_det), 100.0)
        assert got is not None, f"{n_det} did not return within {TIMEOUT_S}s"
        rel = abs(got[0] - ref[0]) / ref[0]
        assert rel < 0.01, (f"{n_det} total {got[0]:.6g} differs from "
                            f"{ref[0]:.6g} by {rel:.3%} -- more than the one "
                            f"pixel of extra detector can account for")


def test_centred_object_is_unaffected():
    """A centred object never trips the bug and must be untouched by the
    fix -- the control that catches a 'fix' which changes physics."""
    got = _run(np.array([1239, 1179]), 0.0)
    assert got is not None, "centred case did not return"
    assert got[1] and got[0] > 0.0


if __name__ == "__main__":
    for fn in (test_missed_rays_terminate_with_odd_detector,
               test_detector_parity_does_not_change_the_result,
               test_centred_object_is_unaffected):
        print(f"{fn.__name__} ... ", end="", flush=True)
        try:
            fn()
            print("PASS")
        except AssertionError as e:
            print(f"FAIL\n    {e}")
