#%% DEMO 03: Generate sample data and add realistic CT noise to it.
#
# This demo will show how to generate sample data for image reconstruction
#
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
#%%
import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise

geo = tigre.geometry_default(high_resolution=False)

#%% Define angles of projection and load phatom image

# define projection angles (in radians)
angles = np.linspace(0, 2 * np.pi, 50)
# load phatnom image
head = sample_loader.load_head_phantom(geo.nVoxel)

# Simulate forward projection.
# To match with mathematical notation, the projection operation is called Ax
projections = tigre.Ax(head, geo, angles)

# Add realistic noise. Adds photon scattering noise ('Poisson') and
# electronic noise of the detector ('Gaussian').
#
# 'Poisson' is related to the maximum photon count in the detector. 1e5 is
# a standard clinical nuber, reduce it for more noise
# 'Gaussian' is related to possible electronic noise in the detector. mean
# of 0 and std of 10 is common in clinical scenario. Increase std for more
# noise.
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

#%% Plot Projections
tigre.plotproj(projections)
# plot noise
tigre.plotproj(projections - noise_projections)
