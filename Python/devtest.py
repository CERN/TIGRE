from __future__ import print_function

import numpy as np
from tigre.geometry import geometry
import tigre_demo_file.Test_data.data_loader as data_loader
from tigre.Ax import Ax
from tigre.Atb import Atb
from tigre.Algorithms.FDK import FDK
from tigre.Algorithms.SIRT import SIRT
from tigre.Utilities.plotproj import plotproj
from tigre.Utilities.plotproj import ppslice
from tigre.Utilities.plotImg import plotImg
from matplotlib import pyplot as plt
#---------------GEOMETRY---------------------------

geo = geometry()
geo.DSD = 1536                                     # Distance Source Detector      (mm)
geo.DSO = 1000                                     # Distance Source Origin        (mm)
# Detector parameters
geo.nDetector = np.array((512, 512))               # number of pixels              (px)
geo.dDetector = np.array((0.8, 0.8))               # size of each pixel            (mm)
geo.sDetector = geo.nDetector * geo.dDetector    # total size of the detector    (mm)
# Image parameters
geo.nVoxel = np.array((256, 256, 256))             # number of voxels              (vx)
geo.sVoxel = np.array((256, 256, 256))             # total size of the image       (mm)
geo.dVoxel = geo.sVoxel/geo.nVoxel               # size of each voxel            (mm)
geo.mode = 'cone'

#geo.nVoxel=np.array([128,128,128])
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
#cube =data_loader.load_cube((64,64,64))


#---------------------ANGLES--------------------------------------------

angles_1 = np.linspace(0, 2*np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100),dtype=np.float32)*np.array(np.pi/4, dtype=np.float32)
angles_3 = np.zeros((100),dtype=np.float32)
angles= np.vstack((angles_1,angles_2,angles_3)).T


proj_1 = Ax(source_img,geo,angles,'interpolated')
geo.check_geo(angles)
#--------------------Back Projection ----------------
#fbp = Atb(proj_1,geo,angles,'FDK')
fbp = FDK(proj_1,geo,angles)


#---------------PLOTS-------------------------------
#plotproj(proj_1)
#plt.matshow(np.hstack((fbp[32],otherfbp.T[32])))

#plt.matshow(abs(fbp[32]-otherfbp.T[32]))

#plt.show()
#plotImg(np.hstack((source_img,fbp)))
plotImg(fbp)
plt.matshow(fbp[32])
plt.show()