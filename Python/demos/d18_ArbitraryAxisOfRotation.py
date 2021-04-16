## DEMO 18: Arbitrary axis of rotation
#
#
#
# Some modenr CT geometires are starting to be a bit more complex, one of
# the common things being arbitrary axis of rotation i.e. the detector and the
# source can move not in a circular path, but in a "spherical" path.
#
# In TIGRE this has been implemented by defining the rotation with 3
# angles, specifically the ZYZ configuration of Euler angles.
#
#  This demo shows how to use it.
#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# This file is part of the TIGRE Toolbox
# # Copyright (c) 2015, University of Bath and
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
#%%Initialize
import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs

#%% Geometry
geo = tigre.geometry_default(high_resolution=False)

#%% Load data and generate projections
# define angles
angles = np.linspace(0, 2 * np.pi, 100)
## Define angles
numProjs = 100


anglesY = np.linspace(0, 2 * np.pi, numProjs)
anglesZ2 = anglesY
anglesZ1 = np.pi * np.sin(np.linspace(0, 2 * np.pi, numProjs))
angles = np.vstack([anglesZ1, anglesY, anglesZ2]).T

## Get Image

head = sample_loader.load_head_phantom(geo.nVoxel)

## Project

projections = tigre.Ax(head, geo, angles)

tigre.plotproj(projections, anglesY)  # angle information not right in the title
## Reconstruct:

# Note, FDK will not work.

imgSIRT = algs.sirt(projections, geo, angles, 50)
imgCGLS = algs.cgls(projections, geo, angles, 10)

tigre.plotimg(np.concatenate([head, imgSIRT, imgCGLS], axis=1), dim="z")
