##% Demo 6: Algorithms01
#
# In this demo the usage of the FDK is explained
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
# Load thorax phatom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections = tigre.Ax(head, geo, angles)
# add noise
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

# %% Usage of FDK

# the FDK algorithm has been taken and modified from
# 3D Cone beam CT (CBCT) projection backprojection FDK, iterative reconstruction Matlab examples
# https://www.mathworks.com/matlabcentral/fileexchange/35548-3d-cone-beam-ct--cbct--projection-backprojection-fdk--iterative-reconstruction-matlab-examples

# The algorithm takes, as eny of them, 3 mandatory inputs:
# PROJECTIONS: Projection data
# GEOMETRY   : Geometry describing the system
# ANGLES     : Propjection angles
# And has a single optional argument:
# FILTER: filter type applied to the projections. Possible options are
#        'ram_lak' (default)
#        'shepp_logan'
#        'cosine'
#        'hamming'
#        'hann'
# The choice of filter will modify the noise and sopme discreatization
# errors, depending on which is chosen.
#
imgFDK1 = algs.fdk(noise_projections, geo, angles, filter="hann")
imgFDK2 = algs.fdk(noise_projections, geo, angles, filter="ram_lak")

# They look quite the same
tigre.plotimg(np.concatenate([imgFDK1, imgFDK2], axis=1), dim="Z")

# but it can be seen that one has bigger errors in the whole image, while
# hte other just in the boundaries
tigre.plotimg(np.concatenate([abs(head - imgFDK1), abs(head - imgFDK2)], axis=1), dim="Z")
