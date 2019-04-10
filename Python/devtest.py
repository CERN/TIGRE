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

#geo1 = tigre.geometry(mode='cone', high_quality=False, default=True)
geo = tigre.geometry(mode='cone', nVoxel=np.array([256,256,256]), default=True)

niter = 10
nangles = 100
angles = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
#head = np.load('src_img_cubic_256.npy') #data_loader.load_head_phantom(geo.nVoxel)
head = data_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(head,geo,angles)
output = algs.sirt(proj,geo,angles,niter)

# print(proj.shape)
# plt.imshow(proj[50])
# plt.show()