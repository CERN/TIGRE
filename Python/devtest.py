from __future__ import print_function
import os
import sys
import tigre

import numpy as np
from matplotlib import pyplot as plt
from tigre.demos.Test_data import data_loader
from tigre.utilities.Ax import Ax
import tigre.algorithms as algs
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
# ---------------GEOMETRY---------------------------

geo = tigre.geometry(mode='parallel',nVoxel = np.array([64,64,64]))
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)


# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PROJECTION----------------------

proj = Ax(source_img,geo,angles)

# ---------------------RECONSTRUCTION------------------
sart = IterativeReconAlg(proj,geo,angles,niter=10,**dict(blocksize=1))
def print_algs(alg):
    for item in alg.__dict__:
        if item == 'geo':
            print('geo')
        if hasattr(getattr(alg,item),'shape'):
            print(item + ' :' + str(getattr(alg,item).shape))
        else:
            print(item,getattr(alg,item))
print_algs(sart)
from tigre.algorithms.pocs_algorithms import ASD_POCS
asd_pocs = ASD_POCS(proj,geo,angles,niter=10)
print_algs(asd_pocs)
res = sart.run_main_iter()

# ---------------------PLOT----------------------------
