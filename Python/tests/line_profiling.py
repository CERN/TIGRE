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

#----------------PROFILINGIMPORT-------------------

from line_profiler import LineProfiler
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

class lineprofileroveride(IterativeReconAlg):
    def __init__(self,proj,geo,angles,niter,**kwargs):
        lp = LineProfiler()
        lp_wrapper = lp(super(lineprofileroveride,self).__init__)
        lp_wrapper(proj,geo,angles,niter,**kwargs)
        lp.print_stats()

alg = lineprofileroveride(proj,geo,angles,10,**dict(blocksize = 10))




#sart = IterativeReconAlg(proj,geo,angles,niter=10, **dict(blocksize =20))
#lp_wrapper = lp(getattr(sart,sart.dataminimizing))
#setattr(sart, sart.dataminimizing, lp_wrapper)
#sart.run_main_iter()
#lp.print_stats()
from tigre.algorithms.pocs_algorithms import ASD_POCS
lp = LineProfiler()
asd_pocs = ASD_POCS(proj,geo,angles,2, **dict(blocksize=20))
lp_wrapper_data_minimizing = lp(getattr(asd_pocs,asd_pocs.dataminimizing))
setattr(asd_pocs, asd_pocs.dataminimizing, lp_wrapper_data_minimizing)
lp_wrapper_main_iter = lp(asd_pocs.run_main_iter)
lp_wrapper_main_iter()
lp.print_stats()

# ---------------------PLOT----------------------------
