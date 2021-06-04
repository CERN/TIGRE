#%% Demo 13: Helical Geometry tests
#
#
# This demo shows an example of TIGRE working on Helical scan geometries
#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# This file is part of the TIGRE Toolbox
#
# Copyright (c) 2015, University of Bath and
#                     CERN-European Organization for Nuclear Research
#                     All rights reserved.
#
# License:            Open Source under BSD.
#                     See the full license at
#                     https://github.com/CERN/TIGRE/blob/master/LICENSE
#
# Contact:            tigre.toolbox@gmail.com
# Codes:              https://github.com/CERN/TIGRE/
# Coded by:           Ander Biguri
# --------------------------------------------------------------------------
##
#%%Initialize
import tigre
import numpy as np
from tigre.utilities import sample_loader
import tigre.algorithms as algs
import time

#%% Geometry
geo = tigre.geometry_default(high_resolution=False)

angles = np.linspace(0, 2 * np.pi, 100)
angles = np.hstack([angles, angles, angles])  # loop 3 times

# Load thorax phatom data
head = sample_loader.load_head_phantom(geo.nVoxel)

# This makes it helical
geo.offOrigin = np.zeros((angles.shape[0], 3))
geo.offOrigin[:, 2] = np.linspace(
    -1024 / 2 + 128, 1024 / 2 - 128, angles.shape[0]
)  # about 256^3 images fit int he detector with this size.


# project data
data = tigre.Ax(head, geo, angles)

# Uncomment if you want to see the data
# plotProj(data,angles);
## Reconstruct Helical

SIRTimg = algs.sirt(data, geo, angles, 30)
# SARTimg=SART(data,geo,angles,30); # takes time
CGLSimg = algs.cgls(data, geo, angles, 20)
## Plot results

# CGLS and SIRT
tigre.plotImg(np.concatenate([head, SIRTimg, CGLSimg], axis=1), dim="z", step=3, clims=[0, 1])
