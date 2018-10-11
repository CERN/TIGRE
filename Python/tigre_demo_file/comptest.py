from __future__ import print_function
import os
import sys
import scipy as sp
import numpy as np
import tigre.geometry as geometry
import Test_data.data_loader as data_loader
from matplotlib import pyplot as plt
import scipy.io
from tigre.Utilities.plotproj import ppslice
from tigre.Utilities.plotproj import plotproj
from tigre.Utilities.plotImg import plotImg
from tigre.Utilities.Measure_Quality import Measure_Quality as MQ
from _Ax import Ax
from _Atb import Atb
from tigre.Algorithms.SART import SART
from tigre.Algorithms.SIRT import SIRT
from tigre.Algorithms.OS_SART import OS_SART
from tigre.Algorithms.FDK import FDK
import time
geo = geometry.TIGREParameters(high_quality=False)

source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
angles = np.linspace(0, 2*np.pi, 100, dtype=np.float32)

projsirt = Ax(source_img, geo, angles, 'ray-voxel')
fdk=FDK(projsirt,geo,angles)
# blocksize=input('blocksize:')
niter = 5
sart=SART(projsirt,geo,angles,niter)

plotImg(np.hstack((fdk,sart)))


