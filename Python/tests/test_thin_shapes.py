"""Thin problems (fan beam, few slices) give the same results as the same rays/voxels of a full 3D problem.

The forward projectors choose their thread-block shape from the detector (a thin detector, e.g. one row for fan
beam, gets blocks along its long side), and the cone backprojectors their voxel column length from the image (a thin
image gets shorter columns). Detectors of 9+ rows and images of 8+ slices keep the original shapes, so each test
compares a thin problem with the same rays or voxels taken from a full-size one computed the original way.
On a thin detector the cone backprojectors also upload all the projections into one texture, instead of chunks of
32: the fan-beam backprojection test compares that with a 17-row detector, uploaded in chunks.
"""
import numpy as np
import pytest

import tigre

pytestmark = pytest.mark.gpu

ANGLES = np.linspace(0, 2 * np.pi, 30, endpoint=False)
# The thin and full problems compute the same rays/voxels from differently rounded positions (and the linear texture
# filtering has 8-bit weights), so they agree to ~1e-6, not bit for bit. A wrong block shape or column would show
# as missing or garbage pixels, orders of magnitude above this.
TOL = 1e-4


def _cone(nvox, ndet, dvox=1.0, ddet=1.5):
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array(nvox))
    geo.DSO, geo.DSD = 300.0, 450.0
    geo.dVoxel = np.array([dvox] * 3, dtype=np.float64)
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.nDetector = np.array(ndet)
    geo.dDetector = np.array([ddet, ddet])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.zeros(3)
    geo.offDetector = np.zeros(2)
    geo.rotDetector = np.zeros(3)
    return geo


def _parallel(nvox, ndet):
    geo = tigre.geometry(mode="parallel", nVoxel=np.array(nvox))
    geo.DSO, geo.DSD = 300.0, 450.0
    geo.dVoxel = np.ones(3)
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.nDetector = np.array(ndet)
    geo.dDetector = np.ones(2)
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.zeros(3)
    geo.offDetector = np.zeros(2)
    return geo


def _image(shape, seed=0):
    rng = np.random.default_rng(seed)
    nz, ny, nx = shape
    z, y, x = np.meshgrid(*[np.linspace(-1, 1, n) if n > 1 else np.zeros(1) for n in shape], indexing="ij")
    img = np.zeros(shape, np.float32)
    for _ in range(6):  # smooth blobs, no symmetry
        c = rng.uniform(-0.5, 0.5, 3)
        w = rng.uniform(0.15, 0.4)
        img += np.exp(-((z - c[0]) ** 2 + (y - c[1]) ** 2 + (x - c[2]) ** 2) / (2 * w * w)).astype(np.float32)
    return img


def _rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-30))


def _matches_one_of(thin, candidates, tol, what):
    errs = [_rel(thin, c) for c in candidates]
    assert min(errs) < tol, "{}: relative difference {} (tolerance {})".format(what, errs, tol)


@pytest.mark.parametrize("ptype", ["Siddon", "interpolated"])
def test_fan_beam_equals_central_row(ptype):
    """A one-slice image on a one-row detector (fan beam) = the central row of a 17-row detector."""
    img = _image((1, 64, 64))
    full = tigre.Ax(img, _cone((1, 64, 64), (17, 96)), ANGLES, ptype)
    thin = tigre.Ax(img, _cone((1, 64, 64), (1, 96)), ANGLES, ptype)
    assert thin.shape == (len(ANGLES), 1, 96)
    _matches_one_of(thin[:, 0, :], [full[:, 8, :]], TOL, "fan beam " + ptype)


@pytest.mark.parametrize("ptype", ["Siddon", "interpolated"])
@pytest.mark.parametrize("rows", [1, 3])
def test_thin_cone_detector_equals_rows(ptype, rows):
    """A detector of 1 or 3 rows, shifted onto rows around j of a 17-row one, gives those rows (3D image).

    Which end the rows are counted from depends on TIGRE's sign of offDetector: both are accepted."""
    img = _image((24, 40, 40))
    nv, dv, j = 17, 1.5, 4  # rows j-rows//2 .. j+rows//2 of the full detector
    full = tigre.Ax(img, _cone((24, 40, 40), (nv, 64)), ANGLES, ptype)
    geo = _cone((24, 40, 40), (rows, 64))
    geo.offDetector = np.array([(j - (nv - 1) / 2) * dv, 0.0])
    thin = tigre.Ax(img, geo, ANGLES, ptype)
    lo = j - rows // 2
    mirror = nv - 1 - j - rows // 2
    _matches_one_of(thin, [full[:, lo:lo + rows, :], full[:, mirror:mirror + rows, :]], TOL,
                    "cone, {} rows, {}".format(rows, ptype))


@pytest.mark.parametrize("ptype", ["Siddon", "interpolated"])
def test_parallel_one_slice_equals_slab_rows(ptype):
    """Parallel beam: a one-slice image on a one-row detector = every row of a stack of identical slices."""
    sl = _image((1, 48, 48))
    stack = np.repeat(sl, 12, axis=0)
    full = tigre.Ax(stack, _parallel((12, 48, 48), (12, 72)), ANGLES, ptype)
    thin = tigre.Ax(sl, _parallel((1, 48, 48), (1, 72)), ANGLES, ptype)
    for r in (0, 5, 11):
        _matches_one_of(thin[:, 0, :], [full[:, r, :]], TOL, "parallel row {} {}".format(r, ptype))


@pytest.mark.parametrize("bptype", ["FDK", "matched"])
@pytest.mark.parametrize("nz", [1, 2, 3, 5])
def test_thin_image_backprojection_equals_slices(bptype, nz):
    """An image of 1-5 slices, placed at slices k.. of a 16-slice one, gets the same backprojection there."""
    n, dz, k = 16, 1.0, 3
    # smooth projections: linear texture interpolation has 8-bit weights, so a sub-1e-5 shift of a sample point can
    # change a weight by 1/256; on smooth data that is far below the tolerance, on noise it would not be
    proj = tigre.Ax(_image((n, 32, 32), seed=2), _cone((n, 32, 32), (24, 48)), ANGLES, "interpolated")
    full = tigre.Atb(proj, _cone((n, 32, 32), (24, 48)), ANGLES, bptype)
    geo = _cone((nz, 32, 32), (24, 48))
    geo.offOrigin = np.array([(k + (nz - 1) / 2 - (n - 1) / 2) * dz, 0.0, 0.0])
    thin = tigre.Atb(proj, geo, ANGLES, bptype)
    assert thin.shape == (nz, 32, 32)
    other = n - k - nz  # the same slices counted from the other end, if offOrigin's z runs the other way
    _matches_one_of(thin, [full[k:k + nz], full[other:other + nz]], TOL, "{} slices, {}".format(nz, bptype))


@pytest.mark.parametrize("bptype", ["FDK", "matched"])
@pytest.mark.parametrize("nangles", [20, 100, 360])
def test_fan_beam_backprojection_equals_tall_detector(bptype, nangles):
    """Fan beam (a one-row detector: all projections in one texture upload) backprojected into one slice = a
    17-row detector whose rows all hold the same values (chunks of 32 projections), into the same slice."""
    angles = np.linspace(0, 2 * np.pi, nangles, endpoint=False)
    fan = _cone((1, 64, 64), (1, 96))
    proj = tigre.Ax(_image((1, 64, 64), seed=3), fan, angles, "interpolated")  # smooth along the detector
    thin = tigre.Atb(proj, fan, angles, bptype)
    full = tigre.Atb(np.ascontiguousarray(np.repeat(proj, 17, axis=1)), _cone((1, 64, 64), (17, 96)), angles, bptype)
    assert thin.shape == (1, 64, 64)
    _matches_one_of(thin, [full], TOL, "fan beam backprojection, {} angles, {}".format(nangles, bptype))
