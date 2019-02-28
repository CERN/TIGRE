from __future__ import print_function
import tigre
import copy
import tigre.algorithms as algs
import numpy as np
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt
import warnings
#warnings.filterwarnings("error")
nVoxel = np.array([64,64,64])
nangles = 100
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.ones((nangles), dtype=np.float32)*np.pi
angles = np.vstack((angles_1, angles_3, angles_3)).T

geo1 = tigre.geometry(nVoxel=nVoxel,mode='parallel',default_geo=True)
geo2 = tigre.geometry(nVoxel=nVoxel,mode='cone',default_geo=True)
geo3 = copy.deepcopy(geo1)
geo3.accuracy = 0.712

print(str(geo1 == geo2))
print(geo1.accuracy)
print(geo3.accuracy)
print(geo1.accuracy ==geo3.accuracy)
print(geo1 == geo1)
print(geo1 == geo3)
"""
src_img = data_loader.load_head_phantom(nVoxel)
proj=tigre.Ax(src_img,geo,angles)
res = tigre.algorithms.fbp(proj,geo,angles)
plt.imshow(res[nVoxel[0]/2])
plt.show()
"""