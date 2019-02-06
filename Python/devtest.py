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

t = time.time()
for filename in os.listdir(os.curdir):
    if filename.endswith('.npy'):
        rm_files += ' ' + filename
        print(rm_files)
os.system('rm' + rm_files)
nVoxel = np.array([64, 64, 64])
nVoxel = np.array([256, 256, 256])

"""
# ---------------PARALLEL GEOMETRY---------------------------

geo_par = tigre.geometry(mode='parallel', nVoxel=nVoxel)
source_img = data_loader.load_head_phantom(number_of_voxels=geo_par.nVoxel)

# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PARALLEL PROJECTION----------------------

proj_par = Ax(source_img, geo_par, angles)

# ---------------------PARALLEL RECONSTRUCTION------------------


alglist = [  # 'sart',
    'sirt',
    #'ossart',
    # 'iterativereconalg',
    #'asd_pocs',
    'fbp',
    # 'cgls'
]

do_algs(alglist, proj_par, geo_par, angles, mode='parallel', niter=20, **dict(nVoxel=nVoxel, blocksize=20))

# ---------------CONE GEOMETRY---------------------------

geo_con = tigre.geometry_default(nVoxel=nVoxel)
source_img = data_loader.load_head_phantom(number_of_voxels=geo_con.nVoxel)

# ---------------------ANGLES-------------------------
nangles = 100
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.ones((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((nangles), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------CONE PROJECTION----------------------
source_img = data_loader.load_head_phantom(number_of_voxels=geo_con.nVoxel)
proj_con = Ax(source_img, geo_con, angles)

# ---------------------CONE RECONSTRUCTION------------------

alglist = [  # 'sart',
    'sirt',
    #'ossart',
    # 'iterativereconalg',
    #'asd_pocs',
    #'FDK',
    # 'cgls'
]

do_algs(alglist, proj_con, geo_con, angles, mode='cone', niter=20, **dict(nVoxel = nVoxel, blocksize=20))

# --------------------------- CGLS for both modes---------------------------------------
do_algs(['cgls'], proj_par, geo_par, angles, mode='Parallel', **dict(nVoxel = [64 ,64 ,64]))
do_algs(['cgls'], proj_con, geo_con, angles, mode='Cone', **dict(nVoxel = [64 ,64 ,64]))
"""
from matplotlib import pyplot as plt

import tigre.algorithms as algs
nangles = 50
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((nangles), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T
geo = tigre.geometry(mode='cone',nVoxel=nVoxel,default_geo=True)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
elapsed = time.time() - t
geo.nDetector = np.array([256, 256])
geo.dDetector = geo.sDetector/geo.nDetector
print(elapsed)

#proj = tigre.Ax(source_img,geo,angles)
from _minTV import minTV

elapsed=[]
for i in range(1):
    t = time.time()
    proj = tigre.Ax(source_img,geo,angles)
    elapsed.append(time.time() - t)


print(np.mean(np.array(elapsed)))
#plt.imshow(proj[5])
#plt.show()
res=tigre.Atb(proj, geo, angles)

#res = algs.awasd_pocs(proj,geo,angles,niter=1,**dict(blocksize=nangles/5))
plt.imshow(res[100])
plt.show()
"""
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
"""
