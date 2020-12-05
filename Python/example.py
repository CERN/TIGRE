from __future__ import print_function
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import tigre
import tigre.algorithms as algs
from tigre.demos.Test_data import data_loader
from tigre.utilities.Measure_Quality import Measure_Quality

# Geometry
#geo1 = tigre.geometry(mode='cone', high_quality=False, default=True)
geo = tigre.geometry(mode='cone', nVoxel=np.array([256,256,256]), default=True)
geo.dDetector = np.array([0.8, 0.8])*2               # size of each pixel            (mm)
geo.sDetector = geo.dDetector * geo.nDetector
# print(geo)

nangles = 100
angles = np.linspace(0, 2 * np.pi, nangles, endpoint=False, dtype=np.float32)

# Prepare projection data
#head = np.load('src_img_cubic_256.npy')
head = data_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(head,geo,angles)

# Reconstruct
niter = 20
fdkout = algs.fdk(proj,geo,angles)
sirtout = algs.ossart(proj,geo,angles,niter,blocksize=20)

# Measure Quality
# 'RMSE', 'MSSIM', 'SSD', 'UQI'
print('RMSE fdk:')
print(Measure_Quality(fdkout,head,['nRMSE']))
print('RMSE ossart')
print(Measure_Quality(sirtout,head,['nRMSE']))

# Plot
fig, axes=plt.subplots(3, 2)
axes[0, 0].set_title('FDK')
axes[0, 0].imshow(fdkout[geo.nVoxel[0]//2])
axes[1, 0].imshow(fdkout[:, geo.nVoxel[1]//2, :])
axes[2, 0].imshow(fdkout[:, :, geo.nVoxel[2]//2])
axes[0, 1].set_title('OS-SART')
axes[0, 1].imshow(sirtout[geo.nVoxel[0]//2])
axes[1, 1].imshow(sirtout[:, geo.nVoxel[1]//2, :])
axes[2, 1].imshow(sirtout[:, :, geo.nVoxel[2]//2])
plt.show()
# tigre.plotProj(proj)
# tigre.plotImg(fdkout)
