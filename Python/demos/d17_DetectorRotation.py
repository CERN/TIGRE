#% DEMO 17: Detector Rotation
#
#
#
#  Some systems have a slight rotation of the detector due to mechanicall
#  inacuracies.
#
#  According to the article "A geometric calibration method for cone beam
#  CT systems" (DOI: 10.1118/1.2198187), only Roll needs to be corrected
#  for in the algorithmic part, as the other 2 can be easily ignored if
#  "sufficiently small". In TIGRE we decided to leave that to the users
#  discretion and implemented the 3 possible detector rotation, per angle
#  if needed.
#
#  This demo shows how to use it.
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

roll = angles
pitch = 0.7 * np.linspace(0, 1, len(angles))
yaw = 0.7 * np.linspace(0, 1, len(angles))

geo.rotDetector = np.array(np.vstack([roll, pitch, yaw])).T

#%%
# Load thorax phatom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections = tigre.Ax(head, geo, angles)
# add noise
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

tigre.plotproj(projections, angles)


#%% recon
imgRotDet = algs.ossart(projections, geo, angles, 50)

#% No rotation
geo.rotDetector = np.array([0, 0, 0])
projections2 = tigre.Ax(head, geo, angles)

imgnoRot = algs.ossart(projections2, geo, angles, 50)
# %% Plot to show that indeed the reconstruction is right

tigre.plotimg(np.concatenate([imgRotDet, imgnoRot], axis=1), dim="z")
