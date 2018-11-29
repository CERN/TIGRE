from __future__ import print_function
import tigre
import numpy as np
import tigre.demos.Test_data.data_loader as data_loader
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from matplotlib import pyplot as plt
import tigre
from tigre.Ax import Ax

# ---------------GEOMETRY---------------------------

geo = tigre.geometry_default(high_quality=False)
source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
geo.mode = 'cone'
print(source_img.shape)
# ---------------------ANGLES-------------------------

angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
angles_3 = np.zeros((100), dtype=np.float32)
angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PROJECTION----------------------
proj_1 = Ax(source_img, geo, angles, 'interpolated')

# --------------------BACK PROJECTION ----------------
alg = IterativeReconAlg(proj_1,geo,angles,2,**dict(verbose=False))
alg.run_main_iter()
sart = alg.getres()
#sart = tigre.algorithms.iterativereconalg(proj_1, geo, angles, 5)
#sart = tigre.algorithms.iterativereconalg(proj_1,geo,angles,niter=5, **dict(blocksize= 1))


#fdk = FDK(proj_1,geo,angles)

# ---------------PLOTS-------------------------------
plt.matshow(sart[32])
plt.show()

