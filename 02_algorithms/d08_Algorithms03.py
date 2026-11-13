"""
Krylov subspace methods
=======================

CGLS, LSQR, LSMR, hybrid LSQR and the GMRES variants.

Source: `Python/demos/d08_Algorithms03.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d08_Algorithms03.py>`__ · `MATLAB/Demos/d08_Algorithms03.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d08_Algorithms03.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d08_Algorithms03.m
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
# Usage CGLS
# ----------

#
#
#  CGLS has the common 4 inputs for iterative algorithms in TIGRE:
#
#  Projections, geometry, angles, and number of iterations
#
# Additionally it contains optional initialization techniques, but we
# recommend not using them. CGLS is already quite fast and using them may
# lead to divergence.
# The options are:
#  'init'    Describes different initialization techniques.
#             -  'none'     : Initializes the image to zeros (default)
#             -  'FDK'      : initializes image to FDK reconstruction
#             -  'multigrid': Initializes image by solving the problem in
#                            small scale and increasing it when relative
#                            convergence is reached.
#             -  'image'    : Initialization using a user specified
#                            image. Not recommended unless you really
#                            know what you are doing.
#  'initimg'    an image for the 'image' initialization. Avoid.

# # use CGLS
imgCGLS, normL2CGLS = algs.cgls(noise_projections, geo, angles, 30, computel2=True)
# use LSQR
imgLSQR, normL2LSQR = algs.lsqr(noise_projections, geo, angles, 30, computel2=True)
# use LSMR
imgLSMR, normL2LSMR = algs.lsmr(noise_projections, geo, angles, 30, computel2=True,lmbda=0)
imgLSMR2, normL2LSMR2 = algs.lsmr(noise_projections, geo, angles, 30, computel2=True,lmbda=30)
# use hybrid LSQR
imghLSQR, normhL2LSQR = algs.hybrid_lsqr(noise_projections, geo, angles, 30, computel2=True)

# AB/BA-GMRES
imgabgmres, normhabgmres = algs.ab_gmres(noise_projections, geo, angles, 30, computel2=True)
imgbagmres, normhbagmres = algs.ba_gmres(noise_projections, geo, angles, 30, computel2=True)
# # AB/BA-GMRES with FDK backprojector
imgabgmresfdk, normhabgmresfdk = algs.ab_gmres(noise_projections, geo, angles, 30, computel2=True,backprojector="FDK")
imgbagmresfdk, normhbagmresfdk = algs.ba_gmres(noise_projections, geo, angles, 30, computel2=True,backprojector="FDK")

# SIRT for comparison.
imgSIRT, normL2SIRT = algs.sirt(noise_projections, geo, angles, 60, computel2=True)

# %%
# plot results
# ------------

#
# We can see that CGLS gets to the same L2 error in less amount of
# iterations.

plt.plot(np.vstack((normL2CGLS[0, :], normL2SIRT[0, 0:30],normL2LSMR[0, :],normL2LSMR2[0, :],normhL2LSQR[0, :],normhabgmres[0,:],normhbagmres[0,:],normhabgmresfdk[0,:],normhbagmresfdk[0,:])).T)
plt.title("L2 error")
plt.xlabel("Iteration")
plt.ylabel("$ |Ax-b| $")
plt.gca().legend(("CGLS","SIRT","LSMR lambda=0","LSMR lambda=30","hybrid LSQR","AB-GMRES","BA-GMRES","AB-GMRES FDK","BA-GMRES FDK"))
plt.show()
# plot images
tigre.plotimg(np.concatenate([np.concatenate([imgCGLS, imgSIRT, imgLSQR,imgabgmres,imgabgmresfdk],axis=1),np.concatenate([imgLSMR, imgLSMR2, imghLSQR, imgbagmres, imgbagmresfdk], axis=1)], axis=2), dim="z", step=2, clims=[0, 2])
# plot errors
tigre.plotimg(np.concatenate([np.concatenate([head-imgCGLS, head-imgSIRT, head-imgLSQR, head-imgabgmres, head-imgabgmresfdk],axis=1),np.concatenate([head-imgLSMR, head-imgLSMR2, head-imghLSQR, head-imgbagmres, head-imgbagmresfdk], axis=1)], axis=2), dim="z", slice=32)
