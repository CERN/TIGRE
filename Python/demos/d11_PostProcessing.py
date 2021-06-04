#%% Demo 11: Postprocessing
#
# This demo demonstrates the available postprocessing tools in TIGRE.by calling the "Measure_Quality.m" function with detailed description.
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
# Coded by:           Manasavee Lohvithee
# --------------------------------------------------------------------------
#%%Initialize
import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs
from matplotlib import pyplot as plt


#%% Geometry
geo = tigre.geometry_default(high_resolution=False)

#%% Load data and generate projections
# define angles
angles = np.linspace(0, 2 * np.pi, 100)
# Load thorax phatom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections = tigre.Ax(head, geo, angles)
# add noise
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

#%% Some recon, FDK for example
imgFDK = algs.fdk(projections, geo, angles)

# TODO, these are not implemented/accesible in python TIGRE
# Issues #270 #271
