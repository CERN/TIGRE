#% DEMO 8:  Algorithms 03. Krylov subspace
#
#
# In this demo the usage of the the Krylov subspace family is explained.
# This family of algorithms iterates trhough the eigenvectors of the
# residual (Ax-b) of the problem in descending order, achieving increased
# convergence rates comparing to SART family.
#
# In cases where the data is good quality, SART type families tend to reach
# to a better image, but when the data gets very big, or has bad quality,
# CGLS is a good and fast algorithm
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

#%% Usage CGLS
#
#
#  CGLS has the common 4 inputs for iterative algorithms in TIGRE:
#
#  Projections, geometry, angles, and number of iterations
#
# Additionally it contains optional initialization tehcniques, but we
# reccomend not using them. CGLS is already quite fast and using them may
# lead to divergence.
# The options are:
#  'Init'    Describes diferent initialization techniques.
#             -  'none'     : Initializes the image to zeros (default)
#             -  'FDK'      : intializes image to FDK reconstrucition
#             -  'multigrid': Initializes image by solving the problem in
#                            small scale and increasing it when relative
#                            convergence is reached.
#             -  'image'    : Initialization using a user specified
#                            image. Not recomended unless you really
#                            know what you are doing.
#  'InitImg'    an image for the 'image' initialization. Avoid.

# use CGLS
imgCGLS, errL2CGLS = algs.cgls(noise_projections, geo, angles, 60, computel2=True)
# SIRT for comparison.
imgSIRT, errL2SIRT = algs.sirt(noise_projections, geo, angles, 60, computel2=True)

#%% plot results
#
# We can see that CGLS gets to the same L2 error in less amount of
# iterations.


plt.plot(np.vstack((errL2CGLS[0, :], errL2SIRT[0, :])).T)
plt.title("L2 error")
plt.xlabel("Iteration")
plt.ylabel("$ log_{10}(|Ax-b|) $")
plt.gca().legend(("CGLS", "SIRT"))
plt.show()
# plot images
tigre.plotimg(np.concatenate([imgCGLS, imgSIRT], axis=1), dim="z", step=2)
# plot errors
tigre.plotimg(np.abs(np.concatenate([head - imgCGLS, head - imgSIRT], axis=1)), dim="z", slice=32)
