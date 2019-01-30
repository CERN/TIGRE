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

nVoxel = np.array([512, 512, 512])


def do_algs(alglist, mode, niter=10, **kwargs):
    for alg in alglist:
        print('Algorithm ' + str(alg).upper() + ' with mode ' + str(geo.mode) + ' and nVoxel ' + str(geo.nVoxel))
        if alg == 'FDK' or alg == 'fbp':
            np.save(alg + ' ' + mode, getattr(algs, alg)(proj, geo, angles, **kwargs))
        else:
            np.save(alg + ' ' + mode, getattr(algs, alg)(proj, geo, angles, niter, **kwargs))


# ---------------PARALLEL GEOMETRY---------------------------

geo = tigre.geometry(mode='parallel', nVoxel=nVoxel)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)

# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PARALLEL PROJECTION----------------------

proj = Ax(source_img, geo, angles)

# ---------------------PARALLEL RECONSTRUCTION------------------

alglist = [ #'sart',
           'sirt',
           'ossart',
           'iterativereconalg',
           'asd_pocs',
           'fbp',
           'cgls']

print(geo.nVoxel)
do_algs(alglist,mode='parallel',niter=20)

# ---------------CONE GEOMETRY---------------------------

geo = tigre.geometry_default(high_quality=True)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)

# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------CONE PROJECTION----------------------

proj = Ax(source_img, geo, angles)

# ---------------------CONE RECONSTRUCTION------------------

alglist = [#'sart',
           'sirt',
           'ossart',
           'iterativereconalg',
           'asd_pocs',
           'FDK',
           'cgls']

do_algs(alglist,mode='cone',niter=20)