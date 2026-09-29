"""Demo 26: Wang detector-offset weights on a short scan.

Wang weighting compensates a displaced detector on a FULL circular scan: it
ramps one side of the detector down and relies on the opposing (beta + pi)
views to bring the coverage back to uniform. On a short scan those opposing
views do not exist, so the ramp survives into the reconstruction as a
one-sided shading (a "teardrop"). Any non-zero offDetector switches the
weights on - a sub-pixel calibrated offset is enough - so a short scan with a
calibrated geometry was silently shaded.

FDK now skips the Wang weights when the angles do not cover a full circle
(the same test MATLAB's FDK uses to decide Parker weighting). This demo
reconstructs a uniform cylinder from a 200-degree scan three ways and reports
the left/right interior ratio:

  centred detector            -> balanced (reference)
  0.1 mm offset, Wang forced  -> the teardrop (what happened before)
  0.1 mm offset, default      -> Wang skipped, balanced again

Wang weights are still applied on full scans, where they belong.
"""
import numpy as np
import tigre
import tigre.algorithms as algs
import matplotlib.pyplot as plt

geo = tigre.geometry_default(high_resolution=False)
geo.nVoxel = np.array([64, 64, 64])
geo.sVoxel = np.array([64.0, 64.0, 64.0])
geo.dVoxel = geo.sVoxel / geo.nVoxel

zz, yy, xx = np.mgrid[:64, :64, :64]
mu = 0.02
phantom = (((yy - 31.5) ** 2 + (xx - 31.5) ** 2 <= 20 ** 2) * mu).astype(np.float32)
core = (yy - 31.5) ** 2 + (xx - 31.5) ** 2 <= 14 ** 2

angles = np.linspace(0, np.deg2rad(200), 200)      # short scan
proj = tigre.Ax(phantom, geo, angles)


def lr(vol):
    left = float(vol[16:48, :, :32][core[16:48, :, :32]].mean())
    right = float(vol[16:48, :, 32:][core[16:48, :, 32:]].mean())
    return left, right


vol_ref = algs.fdk(proj, geo, angles)                 # centred detector

geo_off = tigre.geometry_default(high_resolution=False)
geo_off.nVoxel, geo_off.sVoxel, geo_off.dVoxel = geo.nVoxel, geo.sVoxel, geo.dVoxel
geo_off.offDetector = np.array([0.0, 0.1])            # 0.1 mm: a fraction of a pixel

# The old behaviour: Wang weights applied because offDetector != 0. Reproduce
# it by calling the Wang path directly on a full-circle-shaped call.
from tigre.algorithms import single_pass_algorithms as spa
_orig = spa.is_short_scan
spa.is_short_scan = lambda a: False                   # force the pre-fix path
vol_wang = algs.fdk(proj, geo_off, angles)
spa.is_short_scan = _orig
vol_gated = algs.fdk(proj, geo_off, angles)           # new default: skipped

fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
for ax, vol, title in zip(
        axes, (vol_ref, vol_wang, vol_gated),
        ("centred detector (reference)",
         "0.1 mm offset, Wang weights applied (before)",
         "0.1 mm offset, Wang skipped on the short scan (after)")):
    l, r = lr(vol)
    ax.imshow(vol[32], cmap="gray", vmin=0, vmax=1.5 * mu)
    ax.set_title(f"{title}\nleft/right interior ratio {l / r:.2f}", fontsize=9)
    ax.axis("off")
    print(f"{title}: left {l:.5f} right {r:.5f} ratio {l / r:.3f}")
fig.suptitle("200-degree scan of a uniform cylinder - Wang detector-offset weights need a full circle",
             fontsize=10)
plt.tight_layout()
plt.savefig("d26_short_scan_wang.png", dpi=120)
plt.show()
