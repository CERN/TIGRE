from __future__ import print_function
import tigre
import copy
import tigre.algorithms as algs
import numpy as np
import tigre.algorithms as algs
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt

nVoxel = np.array([256,256,256])
nangles = 100
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.ones((nangles), dtype=np.float32)*np.pi
angles = np.vstack((angles_1, angles_3, angles_3)).T

geo_cone = tigre.geometry(nVoxel=nVoxel,mode='parallel',default_geo=True)
geo_para = tigre.geometry(nVoxel=nVoxel,mode='cone',default_geo=True)

head = data_loader.load_head_phantom(nVoxel)

proj_cone = tigre.Ax(head,geo_cone,angles,'interpolated')
proj_par = tigre.Ax(head,geo_para,angles)


"""
backproj_cone = tigre.Atb(proj_cone,geo_cone,angles,'FDK')
backproj_par = tigre.Atb(proj_par,geo_para,angles,'FDK')
"""
#algossart = algs.ossart(proj_cone,geo_cone,angles,niter=10)

"""
plt.figure()
plt.subplot(2,1,1)
plt.imshow(backproj_cone[128])
plt.subplot(2,1,2)
plt.imshow(backproj_par[128])

plt.show()
"""

