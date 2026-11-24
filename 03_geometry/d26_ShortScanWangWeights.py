"""
Shifted detector on a short scan
================================

Why Wang weights are skipped on short scans, shown on a uniform cylinder.

Source: `Python/demos/d26_ShortScanWangWeights.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d26_ShortScanWangWeights.py>`__ · `MATLAB/Demos/d26_ShortScanWangWeights.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d26_ShortScanWangWeights.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d26_ShortScanWangWeights.m
      :language: matlab
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
