from __future__ import print_function
import tigre
import copy
import tigre.algorithms as algs
import numpy as np
import tigre.algorithms as algs
import time
import sys
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt
from tigre.utilities.Measure_Quality import Measure_Quality
#geo1 = tigre.geometry(mode='cone', high_quality=False, default=True)
geo = tigre.geometry(mode='cone', nVoxel=np.array([256,256,256]),default=True)
geo.dDetector = np.array([0.8, 0.8])*2               # size of each pixel            (mm)
geo.sDetector = geo.dDetector * geo.nDetector

niter = 10
nangles = 100
angles = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
#head = np.load('src_img_cubic_256.npy') #data_loader.load_head_phantom(geo.nVoxel)
head = data_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(head,geo,angles)
fdkout = algs.fdk(proj,geo,angles)
sirtout = algs.ossart(proj,geo,angles,20,blocksize=20)
# 'RMSE'
# 'MSSIM'
# 'SSD'
# 'UQI'
print('RMSE fdk:')
print(Measure_Quality(fdkout,head,['nRMSE']))

print('RMSE ossart')
print(Measure_Quality(sirtout,head,['nRMSE']))
plt.subplot(211)
plt.imshow(fdkout[geo.nVoxel[0]/2])
plt.subplot(212)
plt.imshow(sirtout[geo.nVoxel[0]/2])
plt.colorbar()
plt.show()
# print(proj.shape)
# plt.imshow(proj[50])
# plt.show()