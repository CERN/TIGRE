from __future__ import print_function

import numpy as np
import tigre.geometry_default as geometry
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

geo = geometry.TIGREParameters(high_quality=False)
geo.mode = 'cone'
geo.COR = None
#geo.nVoxel=np.array([128,128,128])
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
#cube =data_loader.load_cube((64,64,64))


#---------------------ANGLES--------------------------------------------

angles_1 = np.linspace(0, 2*np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100),dtype=np.float32)*np.array(np.pi/2, dtype=np.float32)
angles_3 = np.zeros((100),dtype=np.float32)
angles= np.vstack((angles_1,angles_3,angles_3)).T


geo.check_geo(angles)
proj_1 = Ax(source_img,geo,angles,'ray-voxel')

#--------------------Back Projection ----------------
fbp = Atb(proj_1,geo,angles,"FDK")
#fbp = FDK(proj_1,geo,angles)


#---------------PLOTS-------------------------------

print("src_img: "+ str(source_img.shape))
print("nVoxel: "+ str(geo.nVoxel))
plt.matshow(source_img[32])

print("proj: "+ str(proj_1.shape))
plt.matshow(proj_1[50])

print("res: " + str(fbp.shape))
plt.matshow(fbp[32])
#plotImg(fbp)

plt.colorbar()
plt.show()
#plotImg(np.hstack((source_img,fbp)))