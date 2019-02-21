from __future__ import print_function
import tigre
import tigre.algorithms as algs
import numpy as np
from tigre.demos.Test_data import data_loader
from matplotlib import pyplot as plt
import warnings
#warnings.filterwarnings("error")
nVoxel = np.array([32,64,128])
nangles = 100
angles_1 = np.linspace(0, 2 * np.pi, nangles, dtype=np.float32)
angles_2 = np.zeros((nangles), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.ones((nangles), dtype=np.float32)*0.0001
angles = np.vstack((angles_1, angles_3, angles_3)).T
geo = tigre.geometry(mode='parallel',nVoxel=nVoxel,default_geo=True)
geo.accuracy=1.
src_img = data_loader.load_head_phantom(nVoxel)
proj = tigre.Ax(src_img,geo,angles)
algs.cgls(proj,geo,angles,niter=24)
