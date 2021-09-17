##% Demo 7: Algorithms02
#
# In this demo the usage of the algorithms on the SART family among with their options are presented.
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

#%% Algorithms
#
# The main difference between them is the update process.
#   SART: Updates the image projection by projection
#   SIRT: Updates the image using the whole set of projections at once
#   OS-SART: middle ground. Updates the image using a subset of the
#            projections
#
#  Of these algorithms, SART is generally the one reaching a better image
#  (less L2 error) for the same amount of iterations, and SIRT is the
#  worst (still relatively similar). However, SART needs increased
#  computational time per iteration, as it needs to update the image very often,
#  while SIRT only updates the emage ones for the whole sets of projections.
#  OS-SART lies in the middle, reaching similar convergence (L2 error per
#  iteration) than SART but with less computational time than SART.
#
## Usage, with optional parameters.
# In the three algorithms, there are 4 mandatory input arguments:
# Projections, geometry, angles and number of iterations.
#
#
# Optional arguments for all of them
# ==========================================================================
# 'lambda': hyperparameter. The update will be multiplied by this number
# every iteration, to make the steps bigger or smaller. Default: 1
#
lmbda = 1


# 'lambdared': reduction multiplier for the hyperparameter.
# lambda=lambda*lambdared every iterations, so the steps can be smaller
# the further the update. Default=0.99
lambdared = 0.999

# 'Init' : Initialization method. Possible options are
#          None (default). There will be no initialization method, just
#                 the algorithm
#
#          'FDK'  Initialize the image with the result of FDK algorithm
#
#          'multigrid' Initialize using the multigrid method. The image
#                      will be solved in a small scale, and the size of it
#                      will increase until the desired size is reached.
#
#          'image'     Initialzies with a user given image. Not recoomended
#                      unless you really know what you are doing.

initmode = None

# 'InitImg' : related to init. The image to use for initializing the
# algorithm.

# 'verbose': boolean to make the algorithm display (or not) running state.
#            default true.

verbose = True
# 'Quameasopts'  Asks the algorithm for a set of quality measurement
#                parameters. Input should contain a list of desired
#                quality measurement names. Example: ['CC','RMSE','MSSIM']
#                These will be computed in each iteration.
#                Note that this is the change in parameter value per iteration, not the quality of the result.
qualmeas = ["RMSE", "SSD"]

# SIRT and SART both have no extra input parameters.
# =========================================================================
imgSIRT, qualitySIRT = algs.sirt(
    projections,
    geo,
    angles,
    20,
    lmbda=lmbda,
    lmbda_red=lambdared,
    verbose=verbose,
    Quameasopts=qualmeas,
    computel2=True,
)
imgSART, qualitySART = algs.sart(
    projections,
    geo,
    angles,
    20,
    lmbda=lmbda,
    lmbda_red=lambdared,
    verbose=verbose,
    Quameasopts=qualmeas,
    computel2=True,
)

# OS-SART
# ========================================================================
# Additionally OS-SART includes a couple of other parameters, related to
# the subsets.
#
#   'blocksize':   Sets the projection block size used simultaneously. If
#                  BlockSize = 1 OS-SART becomes SART and if  BlockSize = size(angles,2)
#                  then OS-SART becomes SIRT. Default is 20.
blcks = 10
# 'OrderStrategy':  Chooses the subset ordering strategy. Options are
#                  'ordered' :uses them in the input order, but divided
#                  'random'  : orders them randomply
#                  'angularDistance': chooses the next subset with the
#                                     biggest angular distance with the
#                                     ones used.  (default)
order = "random"
imgOSSART, qualityOSSART = algs.ossart(
    projections,
    geo,
    angles,
    20,
    lmbda=lmbda,
    lmbda_red=lambdared,
    verbose=verbose,
    Quameasopts=qualmeas,
    blocksize=blcks,
    OrderStrategy=order,
    computel2=True,
)


#%% Lets have a brief show of the results

plt.subplot(211)
plt.plot(np.log10(np.vstack((qualitySIRT[0, :], qualityOSSART[0, :], qualitySART[0, :]))).T)
plt.title("Convergence")
plt.xlabel("Iteration")
plt.ylabel("$ log_{10}(|Ax-b|) $")
plt.gca().legend(("SIRT", "OS-SART", "SART"))
plt.subplot(212)
plt.plot(np.log10(np.vstack((qualitySIRT[1, :], qualityOSSART[1, :], qualitySART[1, :]))).T)
plt.title("Evolution of RMSE")
plt.gca().legend(("SIRT", "OS-SART", "SART"))
plt.xlabel("Iteration")
plt.ylabel("$ log_{10}(RMSE) $")
plt.show()

#%% plot the results

# It is clear that SART will get to better results for the same amoutn of
# iterations, however, it takes x7 more time to run.

# SART
# OS-SART
# SIRT

tigre.plotimg(np.concatenate([imgSIRT, imgOSSART, imgSART], axis=1), dim="Z", savegif="sarts.gif")

# plot error
tigre.plotimg(
    np.abs(np.concatenate([head - imgSIRT, head - imgOSSART, head - imgSART], axis=1)), dim="Z"
)
