from __future__ import division
from __future__ import print_function

import numpy as np
import tigre
import tigre.algorithms as algs
from tigre.utilities import sample_loader
from tigre.utilities.Measure_Quality import Measure_Quality
import tigre.utilities.gpu as gpu
import matplotlib.pyplot as plt

### This is just a basic example of very few TIGRE functionallity.
# We hihgly recomend checking the Demos folder, where most if not all features of tigre are demoed.


listGpuNames = gpu.getGpuNames()
if len(listGpuNames) == 0:
    print("Error: No gpu found")
else:
    for id in range(len(listGpuNames)):
        print("{}: {}".format(id, listGpuNames[id]))

gpuids = gpu.getGpuIds(listGpuNames[0])
print(gpuids)

# Geometry
# geo1 = tigre.geometry(mode='cone', high_resolution=False, default=True)
geo = tigre.geometry(mode="cone", nVoxel=np.array([256, 256, 256]), default=True)
geo.dDetector = np.array([0.8, 0.8]) * 2  # size of each pixel            (mm)
geo.sDetector = geo.dDetector * geo.nDetector
# print(geo)

nangles = 100
angles = np.linspace(0, 2 * np.pi, nangles, endpoint=False, dtype=np.float32)

# Prepare projection data
head = sample_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(head, geo, angles, gpuids=gpuids)

# Reconstruct
niter = 20
fdkout = algs.fdk(proj, geo, angles, gpuids=gpuids)
ossart = algs.ossart(proj, geo, angles, niter, blocksize=20, gpuids=gpuids)

# Measure Quality
# 'RMSE', 'MSSIM', 'SSD', 'UQI'
print("RMSE fdk:")
print(Measure_Quality(fdkout, head, ["nRMSE"]))
print("RMSE ossart")
print(Measure_Quality(ossart, head, ["nRMSE"]))

# Plot
fig, axes = plt.subplots(3, 2)
axes[0, 0].set_title("FDK")
axes[0, 0].imshow(fdkout[geo.nVoxel[0] // 2])
axes[1, 0].imshow(fdkout[:, geo.nVoxel[1] // 2, :])
axes[2, 0].imshow(fdkout[:, :, geo.nVoxel[2] // 2])
axes[0, 1].set_title("OS-SART")
axes[0, 1].imshow(ossart[geo.nVoxel[0] // 2])
axes[1, 1].imshow(ossart[:, geo.nVoxel[1] // 2, :])
axes[2, 1].imshow(ossart[:, :, geo.nVoxel[2] // 2])
plt.show()
# tigre.plotProj(proj)
# tigre.plotImg(fdkout)
