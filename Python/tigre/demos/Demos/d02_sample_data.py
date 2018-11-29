import os
import sys
# TODO: this is quite nasty; it would be nice to reorganise file structure later so top level folder is always in path
currDir = os.path.dirname(os.path.realpath(__file__))
rootDir = os.path.abspath(os.path.join(currDir, '..'))
if rootDir not in sys.path:  # add parent dir to paths
    sys.path.append(rootDir)

import numpy as np
import tigre.geometry_default as geometry
import Test_data.data_loader as data_loader
import scipy.io
from tigre.utilities.plotproj import plotproj

from _Ax import Ax

TIGRE_parameters = geometry(high_quality=False)

head = data_loader.load_head_phantom(number_of_voxels=TIGRE_parameters.nVoxel)

angles = np.linspace(0, 2*np.pi, 20, dtype=np.float32)
# angles = np.linspace(0, 0, 1, dtype=np.float32)
# What is expected to happen here is the memory allocated for the head image is given to the C functions to operate on
# Whilst it will return newly allocated PYTHON(!) memory for each projection; expected to be allocated from C and parsed
# back into python

# I believe it should be ok to just parse geometry each time until we find a way of passing around a C array in python

projections = Ax(head, TIGRE_parameters, angles, 'interpolated')

plotproj(projections)

# m = {'projections': projections}
# scipy.io.savemat('projections_python', m)
