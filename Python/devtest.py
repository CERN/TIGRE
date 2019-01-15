from __future__ import print_function
import os
import sys
"""
import tigre
#import numpy as np
#import tigre.demos.Test_data.data_loader as data_loader
from tigre.algorithms.conjugate_gradient_algorithms import CGLS
from tigre.algorithms.pocs_algorithms import ASD_POCS
from tigre.algorithms.art_family_algorithms import SART
from matplotlib import pyplot as plt
#import tigre
import tigre.algorithms as algs
from tigre.utilities.Ax import Ax
"""

# ---------------GEOMETRY---------------------------

#geo = tigre.geometry(mode='parallel',nVoxel = np.array([64,64,64],dtype=np.float32))
#source_img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)


# ---------------------ANGLES-------------------------

#angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
#angles_2 = np.ones((100), dtype=np.float32) * np.array(np.pi / 4, dtype=np.float32)
#angles_3 = np.zeros((100), dtype=np.float32)
#angles = np.vstack((angles_1, angles_3, angles_3)).T

# --------------------PROJECTION----------------------
def include_headers(filename_list, sdist=False):
    """add hpp and h files to list if sdist is called"""
    c_extensions = ['.cu',".c", ".C", ".cc", ".cpp", ".cxx", ".c++"]
    header_list = []
    for filename in filename_list:
        if sdist:
            header = list(os.path.splitext(filename))
            if header[1] in c_extensions:
                header[1] = '.hpp'
                header_list.append(''.join(header))
    if sdist:
        filename_list.append('tigre/Source/types.hpp')
    return filename_list + header_list





test = ['tigre/Source/projection.cpp',
'tigre/Source/Siddon_projection.cu', 'tigre/Source/Siddon_projection_parallel.cu',
'tigre/Source/ray_interpolated_projection.cu', 'tigre/Source/ray_interpolated_projection_parallel.cu',
'tigre/Source/_types.pxd',
'tigre/Source/_Ax.pyx']

print(include_headers(test,sdist=True))