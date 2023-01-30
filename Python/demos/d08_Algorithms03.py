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

# # use CGLS
imgCGLS, normL2CGLS = algs.cgls(noise_projections, geo, angles, 30, computel2=True)
# use LSQR
imgLSQR, normL2LSQR = algs.lsqr(noise_projections, geo, angles, 30, computel2=True)
# use LSMR
imgLSMR, normL2LSMR = algs.lsmr(noise_projections, geo, angles, 30, computel2=True,lmbda=0)
imgLSMR2, normL2LSMR2 = algs.lsmr(noise_projections, geo, angles, 30, computel2=True,lmbda=30)
# use LSQR
imghLSQR, normhL2LSQR = algs.hybrid_lsqr(noise_projections, geo, angles, 30, computel2=True)

# AB/BA-GMRES
imgabgmres, normhabgmres = algs.ab_gmres(noise_projections, geo, angles, 30, computel2=True)
imgbagmres, normhbagmres = algs.ba_gmres(noise_projections, geo, angles, 30, computel2=True)
# # AB/BA-GMRES with FDK backprojector
imgabgmresfdk, normhabgmresfdk = algs.ab_gmres(noise_projections, geo, angles, 30, computel2=True,backprojector="FDK")
imgbagmresfdk, normhbagmresfdk = algs.ba_gmres(noise_projections, geo, angles, 30, computel2=True,backprojector="FDK")


# SIRT for comparison.
imgSIRT, normL2SIRT = algs.sirt(noise_projections, geo, angles, 60, computel2=True)

#%% plot results
#
# We can see that CGLS gets to the same L2 error in less amount of
# iterations.



plt.plot(np.vstack((normL2CGLS[0, :], normL2SIRT[0, 0:30],normL2LSMR[0, :],normL2LSMR2[0, :],normhL2LSQR[0, :],normhabgmres[0,:],normhbagmres[0,:],normhabgmresfdk[0,:],normhbagmresfdk[0,:])).T)
plt.title("L2 error")
plt.xlabel("Iteration")
plt.ylabel("$ |Ax-b| $")
plt.gca().legend(("CGLS", "SIRT","LSMR lambda=0", "LSMR lambda=30","hybrid LSQR","AB-GMRES","BA-GMRES","AB-GMRES FDK","BA-GMRES FDK"))
plt.show()
# plot images
tigre.plotimg(np.concatenate([np.concatenate([imgCGLS, imgSIRT, imgLSQR,imgabgmres,imgabgmresfdk],axis=1),np.concatenate([imgLSMR, imgLSMR2, imghLSQR,imgbagmres,imgbagmresfdk], axis=1)], axis=2), dim="z", step=2,clims=[0, 2])
# plot errors
tigre.plotimg(np.concatenate([np.concatenate([head-imgCGLS, head-imgSIRT, head-imgLSQR, head-imgabgmres, head-imgabgmresfdk],axis=1),np.concatenate([head-imgLSMR, head-imgLSMR2, head-imghLSQR,head-imgbagmres,head-imgbagmresfdk], axis=1)], axis=2), dim="z", slice=32)
