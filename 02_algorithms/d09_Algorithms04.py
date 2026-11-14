"""
Total variation algorithms
==========================

ASD-POCS and its variants, and their parameters.

Source: `Python/demos/d09_Algorithms04.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d09_Algorithms04.py>`__ · `MATLAB/Demos/d09_Algorithms04.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d09_Algorithms04.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d09_Algorithms04.m
      :language: matlab
"""

# %%
# Initialize
# ----------

import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs
from matplotlib import pyplot as plt
from tigre.utilities.im3Dnorm import im3DNORM

# %%
# Geometry
# --------

geo = tigre.geometry_default(high_resolution=False)

# %%
# Load data and generate projections
# ----------------------------------

# define angles
angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
# Load thorax phantom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections = tigre.Ax(head, geo, angles)
# add noise
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

# %%
# Lets create a OS-SART test for comparison
# -----------------------------------------

imgOSSART = algs.ossart(noise_projections, geo, angles, 10)

# %%
# Total Variation algorithms
# --------------------------

#
#  ASD-POCS: Adaptative Steepest Descent-Projection On Convex Subsets
# Often called POCS
# ==========================================================================
# ==========================================================================
#  ASD-POCS minimizes At-B and the TV norm separately in each iteration,
#  i.e. reconstructs the image first and then reduces the TV norm, every
#  iteration. As the other algorithms the mandatory inputs are projections,
#  geometry, angles and maximum iterations.
#
# ASD-POCS has a variety of optional arguments, and some of them are crucial
# to determine the behavior of the algorithm. The advantage of ASD-POCS is
# the power to create good images from bad data, but it needs a lot of
# tuning.
#
#  Optional parameters that are very relevant:
# ----------------------------------------------
#    'maxl2err'    Maximum L2 error to accept an image as valid. This
#                  parameter is crucial for the algorithm, determines at
#                  what point an image should not be updated further.
#                  Default is 20% of the FDK L2 norm.
#
# it's called epsilon in the paper
epsilon = (
    im3DNORM(tigre.Ax(algs.fdk(noise_projections, geo, angles), geo, angles) - noise_projections, 2)
    * 0.15
)
#   'alpha':       Defines the TV hyperparameter. Default is 0.002.
#                  However the paper mentions 0.2 as a good choice
alpha = 0.002

#   'tviter':      Defines the amount of TV iterations performed per SART
#                  iteration. Default is 20

ng = 25

# Other optional parameters
# ----------------------------------------------
#   'lmbda':      Sets the value of the hyperparameter lambda for the SART
#                 iterations. Default is 1
#
#   'lmbda_red':   Reduction of lambda Every iteration
#                  lmbda=lmbda_red*lmbda. Default is 0.99
#
lmbda = 1
lambdared = 0.9999  # you generally want 1

#   'alpha_red':   Defines the reduction rate of the TV hyperparameter
#                  alpha=alpha_red*alpha. Default is 0.95
alpha_red = 0.95

#   'rmax':       The maximum allowed image/TV update ratio. If the TV
#                  update changes the image more than this, the parameter
#                  will be reduced. Default is 0.95
ratio = 0.94

#   'verbose'      1 or 0. Default is 1. Gives information about the
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

imgPICCS = algs.piccs(
    noise_projections,
    geo,
    angles,
    10,
    prior = imgOSSART,
    tviter=ng,
    maxl2err=epsilon,
    alpha=alpha,  # less important.
    lmbda=lmbda,
    lmbda_red=lambdared,
    rmax=ratio,
    verbose=verb,
    prior_ratio = 0.5,
)

#  OS_ASD_POCS: Ordered Subset-TV algorithm
# ==========================================================================
# ==========================================================================
#
# The logical next step to improve ASD-POCS is substituting SART with a
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
# This is a more edge preserving algorithms than ASD_POCS.
# delta is the cut-off value of a normalized edge exponential weight....
# not super clear, but it controls at which point you accept something as real vs noise edge.

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

#  IRN-CGLS-TV: CGLS with TV regularization
# ==========================================================================
# ==========================================================================
#
# Should converge fairly fast to a TV regularized solution
imgCGLSTV=algs.irn_tv_cgls(noise_projections,geo,angles,10,lmbda=5,niter_outer=2)

# %%
# plot results
# ------------

# plot images
tigre.plotimg(
    np.concatenate([imgAWASDPOCS, imgOSASDPOCS, imgASDPOCS, imgOSSART, imgCGLSTV, imgPICCS], axis=1), dim="z", step=2
)
