from __future__ import print_function

import numpy as np
import tigre.geometry_default as geometry
import tigre_demo_file.Test_data.data_loader as data_loader
from tigre.Ax import Ax
from tigre.Algorithms.SIRT import SIRT
from tigre.Utilities.plotImg import plotImg

#---------------GEOMETRY---------------------------

geo = geometry.TIGREParameters(high_quality=False)
geo.mode = 'paralell'
geo.COR = None
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)

#---------------------ANGLES--------------------------------------------

angles_1 = np.linspace(0, 2*np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100),dtype=np.float32)*np.array(np.pi/2, dtype=np.float32)
angles_3 = np.zeros((100),dtype=np.float32)
angles= np.vstack((angles_1,angles_2,angles_3)).T


geo.check_geo(angles)

proj_1 = Ax(source_img,geo,angles,'ray-voxel')

#--------------------Back Projection ----------------
fbp = SIRT(proj_1,geo,angles,10,'FDK')
plotImg(fbp)