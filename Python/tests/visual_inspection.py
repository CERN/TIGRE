from __future__ import print_function
import os
import sys

import tigre
import numpy as np
import tigre.demos.Test_data.data_loader as data_loader
import tigre.algorithms as algs
from tigre.utilities.Ax import Ax
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from matplotlib import pyplot as plt

__doc__ == """
 mode: parallel, cone
 algs: list of algs, seperated by space, to be produced. 
 niter: number of iterations 
 kwargs: list of kwargs to pass to algorithms
"""
# ---------------GEOMETRY---------------------------

geo = tigre.geometry(mode=sys.argv[1],nVoxel = np.array([64,64,64],dtype=np.float32))
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)


# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PROJECTION----------------------

proj = Ax(source_img,geo,angles)

# ---------------------PLOT---------------------------
def plot_algs(alglist, niter = 10, *args):
    l = []
    kwargs = {}
    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
        kwargs.update({k:v})
    for alg in alglist:

            l.append(getattr(algs,alg)(proj,geo,angles,niter,**kwargs)[32])
    res = l[0]
    for mat in l[1:]:
        res= np.vstack((res,mat))
    plt.matshow(res)
    plt.title(alglist)
    if kwargs.has_key('colorbar'):
        plt.colorbar()
    plt.show()



if __name__ == '__main__':
    niter = 10
    args = []
    if len(sys.argv) >= 4:
        niter = int(sys.argv[3])
        args = sys.argv[4:]
    if sys.argv[2] == 'all':
        if sys.argv[1] is 'cone':
            geo = tigre.geometry_default(high_quality=False)
            alglist = ['sart',
                    'sirt',
                    'ossart',
                    'iterativereconalg',
                    'FDK',
                    'asd_pocs',
                    'cgls']
        else:
            alglist = ['sart',
                    'sirt',
                    'ossart',
                    'iterativereconalg',
                    'fbp',
                    'asd_pocs',
                    'cgls']
        plot_algs(alglist, niter, *args)
    else:
        plot_algs(sys.argv[2].split(),niter,*args)