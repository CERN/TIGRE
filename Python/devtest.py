from __future__ import print_function
import tigre

import numpy as np
from tigre.demos.Test_data import data_loader
from tigre.utilities.Ax import Ax
import os
import sys
import time
rm_files = ''
from tests.visual_inspection import plot_algs as do_algs


from matplotlib import pyplot as plt

nVoxel = [64,64,64]
nangles = 50
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((nangles), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T
geo = tigre.geometry(mode='cone',nVoxel=nVoxel,default_geo=True)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
proj = tigre.Ax(source_img,geo,angles)
res = tigre.Atb(proj,geo,angles)
plt.imshow(res[32])
plt.show()

