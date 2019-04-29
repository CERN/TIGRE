from __future__ import print_function
import tigre
import tigre.algorithms as algs
import numpy as np
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt

geo = tigre.geometry(mode='cone', nVoxel=np.array([256, 256, 256]), default_geo=True)
niter = 10
nangles = 100
angles = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)

head = data_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(head,geo,angles)
output = algs.sirt(proj,geo,angles,niter=50)

plt.imshow(output[geo.nVoxel[0]/2])
plt.colorbar()
plt.show()