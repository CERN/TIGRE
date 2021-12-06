#% DEMO 9:   Algorithms 04. Total Variation minimization algorithms
#
#
#  This demo presents the Total variation algorithms in TIGRE. Total
#  variation algorithms try to minimize the variation (gradient) of the
#  image, assuming its piecewise smooth, as most things in nature are (i.e.
#  human body).
#
# This set of algorithms is specially good performing when the noise is
# very big or the number of projections is small, however, they require more
# computational time and memory than the other algorithms to run.
#
#
# This Demo is not fully complete, as several algorithms are not shown
# here, yet are available in TIGRE, namely:
#
# - AwASD_POCS:  Edge preserving ASD_POCS (Adaptative weigthed).
# - OS_AwASD_POCS: OS-version of the previous algorithm
# - PCSD: A version of ASD_POCS that heuristically select some of the
#         parameters, particularly epsilon (maxL2norm)
# - AwPCSD: Edge preserving version of the previous algorithm
# - OS_AwPCSD: block-wise version of the previous
#
# You can use these algorithms in a very similar manner as below examples.
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
from tigre.utilities.im3Dnorm import im3DNORM

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

#%% Lets create a OS-SART test for comparison
imgOSSART = algs.ossart(noise_projections, geo, angles, 10)

#%% Total Variation algorithms
#
#  ASD-POCS: Adaptative Steeppest Descent-Projection On Convex Subsets
# Often called POCS
# ==========================================================================
# ==========================================================================
#  ASD-POCS minimizes At-B and the TV norm separately in each iteration,
#  i.e. reconstructs the image first and then reduces the TV norm, every
#  iteration. As the other algorithms the mandatory inputs are projections,
#  geometry, angles and maximum iterations.
#
# ASD-POCS has a veriety of optional arguments, and some of them are crucial
# to determine the behaviour of the algorithm. The advantage of ASD-POCS is
# the power to create good images from bad data, but it needs a lot of
# tunning.
#
#  Optional parameters that are very relevant:
# ----------------------------------------------
#    'maxL2err'    Maximum L2 error to accept an image as valid. This
#                  parameter is crucial for the algorithm, determines at
#                  what point an image should not be updated further.
#                  Default is 20% of the FDK L2 norm.
#
# its called epsilon in the paper
epsilon = (
    im3DNORM(tigre.Ax(algs.fdk(noise_projections, geo, angles), geo, angles) - noise_projections, 2)
    * 0.15
)
#   'alpha':       Defines the TV hyperparameter. default is 0.002.
#                  However the paper mentions 0.2 as good choice
alpha = 0.002

#   'tviter':      Defines the amount of TV iterations performed per SART
#                  iteration. Default is 20

ng = 25

# Other optional parameters
# ----------------------------------------------
#   'lambda':      Sets the value of the hyperparameter for the SART iterations.
#                  Default is 1
#
#   'lambdared':   Reduction of lambda Every iteration
#                  lambda=lambdared*lambda. Default is 0.99
#
lmbda = 1
lambdared = 0.9999  # you generally want 1


#   'alpha_red':   Defines the reduction rate of the TV hyperparameter
alpha_red = 0.95

#   'Ratio':       The maximum allowed image/TV update ration. If the TV
#                  update changes the image more than this, the parameter
#                  will be reduced.default is 0.95
ratio = 0.94

#   'Verbose'      1 or 0. Default is 1. Gives information about the
#                  progress of the algorithm.

verb = True

imgASDPOCS = algs.asd_pocs(
    noise_projections,
    geo,
    angles,
    10,  # these are very important
    tviter=ng,
    maxl2err=epsilon,
    alpha=alpha,  # less important.
    lmbda=lmbda,
    lmbda_red=lambdared,
    rmax=ratio,
    verbose=verb,
)


#  OS_ASD_POCS: Odered Subset-TV algorithm
# ==========================================================================
# ==========================================================================
#
# The logical next step to imporce ASD-POCS is substituting SART with a
# faster algorithm, such as OS-SART
#
# The parameters are the same as in ASD-POCS, but also have 'BlockSize'


imgOSASDPOCS = algs.os_asd_pocs(
    noise_projections,
    geo,
    angles,
    10,  # these are very important
    tviter=ng,
    maxl2err=epsilon,
    alpha=alpha,  # less important.
    lmbda=lmbda,
    lmbda_red=lambdared,
    rmax=ratio,
    verbose=verb,
    # OSC params
    blocksize=10,
)

#  AwASD_POCS: adaptative weighted ASD_POCS
# ==========================================================================
# ==========================================================================
#
# This is a more edge preserving algorithms than ASD_POCS, in theory.
# delta is the cuttof vlaue of anromalized edge exponential weight....
# not super clear, but it cotnrols at which point you accept something as real vs noise edge.

imgAWASDPOCS = algs.awasd_pocs(
    noise_projections,
    geo,
    angles,
    10,  # these are very important
    tviter=ng,
    maxl2err=epsilon,
    alpha=alpha,  # less important.
    lmbda=lmbda,
    lmbda_red=lambdared,
    rmax=ratio,
    verbose=verb,  # AwASD_POCS params
    delta=np.array([-0.005]),
)

#%% plot results

# plot images
tigre.plotimg(
    np.concatenate([imgAWASDPOCS, imgOSASDPOCS, imgASDPOCS, imgOSSART], axis=1), dim="z", step=2
)
