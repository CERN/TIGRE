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
import os
import sys
rm_files = ''
for filename in os.listdir(os.curdir):
    if filename.endswith('.npy'):
        rm_files += ' ' + filename
        print(rm_files)
os.system('rm' + rm_files)
nVoxel = np.array([64,64,64])


def do_algs(alglist, mode, niter=10, **kwargs):
    for alg in alglist:
        print('Algorithm ' + str(alg).upper() + ' with mode ' + str(geo.mode) + ' and nVoxel ' + str(geo.nVoxel))
        if alg == 'FDK' or alg == 'fbp':
            np.save(alg + '_' + mode, getattr(algs, alg)(proj, geo, angles, **kwargs))
        else:
            np.save(alg + '_' + mode, getattr(algs, alg)(proj, geo, angles, niter, **kwargs))


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

#do_algs(alglist,mode='parallel',niter=20,**dict(blocksize=20))

# ---------------CONE GEOMETRY---------------------------

geo = tigre.geometry_default(high_quality=False)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)

# ---------------------ANGLES-------------------------
nangles=100
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.ones((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((nangles), dtype=np.float32)
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

#do_algs(alglist,mode='cone',niter=20,**dict(blocksize=20))
from tigre.utilities.Atb import Atb
Atb(proj,geo,angles,'FDK')[32]
#plt.colorbar()
#plt.title('FDK')
mat = Atb(proj,geo,angles,'matched')
plt.subplot(221)
plt.imshow(mat[32])
plt.subplot(222)
plt.imshow(mat[:,32])
plt.subplot(223)
plt.imshow(mat[:,:,32])
plt.colorbar()
plt.show()