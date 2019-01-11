from __future__ import print_function
import tigre
import numpy as np
import tigre.demos.Test_data.data_loader as data_loader
from tigre.algorithms.conjugate_gradient_algorithms import CGLS
from tigre.algorithms.pocs_algorithms import ASD_POCS
from tigre.algorithms.art_family_algorithms import SART
from matplotlib import pyplot as plt
import tigre
import tigre.algorithms as algs
from tigre.utilities.Ax import Ax

# ---------------GEOMETRY---------------------------

geo = tigre.geometry(mode='parallel',nVoxel = np.array([64,64,64],dtype=np.float32))
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)


# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PROJECTION----------------------
from tigre.demos.__test import test
test()
"""
# --------------------BACK PROJECTION ----------------
#res = algs.awasd_pocs(proj_1,geo,angles,10,**dict(blocksize=20))
#res = algs.ossart(proj_1,geo,angles,15,**dict(blocksize=12))
res = algs.awasd_pocs(proj_1,geo,angles,5)

# ---------------PLOTS------------------------------

plt.matshow(abs(new_head[32] - source_img[32]))
plt.colorbar()
plt.show()
"""
