from __future__ import print_function
import tigre

import numpy as np
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt
nVoxel = np.array([64,64,64])
nangles = 50
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.ones((nangles), dtype=np.float32)*0.0001
angles = np.vstack((angles_1, angles_3, angles_2)).T
geo = tigre.geometry(mode='parallel',nVoxel=nVoxel,default_geo=True)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
proj = tigre.Ax(source_img,geo,angles)
niter=20
import tigre.algorithms as algs
res = algs.cgls(proj,geo,angles,niter=20,blocksize=12)
plt.imshow(res[32])
plt.colorbar()
plt.show()
#tigre.plotimg(res)
"""
for alg in algs.__all__:
    if alg not in ['fbp','FDK','sart']:
        plt.imshow(getattr(algs,alg)(proj,geo,angles,20,**dict(blocksize=20))[nVoxel[0]/2])
        plt.title(alg)
        plt.colorbar()
        plt.figure()
plt.imshow(algs.FDK(proj,geo,angles)[nVoxel[0]/2])
plt.title('fdk')
plt.colorbar()
plt.show()
"""