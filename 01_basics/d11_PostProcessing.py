"""
Post-process a reconstruction
=============================

Denoise a reconstruction with total variation after it is computed.

Source: `Python/demos/d11_PostProcessing.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d11_PostProcessing.py>`__ · `MATLAB/Demos/d11_PostProcessing.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d11_PostProcessing.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d11_PostProcessing.m
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
from tigre.utilities.im_3d_denoise import im3ddenoise
from tigre.utilities.crop_CBCT import cropCBCT

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
# Some recon, FDK for example
# ---------------------------

imgFDK = algs.fdk(projections, geo, angles)

# %%
# Postprocessing
# --------------

#
#  Currently the postprocessing steps in TIGRE are limited, but there are 2
#  functions you can use
#
# im3ddenoise : Denoises a 3D image, using Total Variation denoising.
#
#  Arguments are the number of iterations and the hyperparameter (same as in
#  SART-TV)
imgdenoised = im3ddenoise(imgFDK, 100, 15)

# cropCBCT: Crops all the sections that lie outide the area that is covered
# by the cone
# Note this does not consider offsets nor does it remove the top and bottom
# cone sections
imcropped = cropCBCT(imgFDK)

# %%
# plot results
# ------------

# denoised image is clearer
tigre.plotImg(
    np.concatenate([imgFDK, imgdenoised, imcropped], axis=1), dim="z"
    )
# however, the denoising has no knoledge of the original data (projections)
# this it doesnt reduce the error. The error increases, specially in small
# areas
tigre.plotImg(
    abs(np.concatenate([head - imgFDK, head - imgdenoised, head - imcropped], axis=1)), dim="z"
    )
