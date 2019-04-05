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

geo = tigre.geometry(mode='cone', nVoxel=np.array([512, 512, 512]), default_geo=True)


niter = 10
nangles = 64
angles = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
#angles_2 = np.zeros((nangles, ), dtype=np.float32)
#angles_3 = np.ones((nangles, ), dtype=np.float32)
#angles = np.vstack((angles_1, angles_2, angles_3)).T
head = data_loader.load_head_phantom(geo.nVoxel)
print(type(head[0][0][0]))


proj = tigre.Ax(head,geo,angles)


#img=algs.sirt(proj,geo,angles,10)

print(sys.getrefcount(proj))
print(proj.sum())
if abs(proj.sum()-41344972)>2:
   print('Unexpected')
plt.figure()
plt.imshow(proj[54])
plt.show()

#tigre.plotproj(proj)
"""
output = tigre.Atb(proj,geo,angles,'FDK')
plt.imshow(output[geo.nVoxel[0]/2])
plt.figure()
plt.imshow(output[:,geo.nVoxel[1]/2])
plt.figure()
plt.imshow(output[:,:,geo.nVoxel[2]/2])
plt.show()
"""
