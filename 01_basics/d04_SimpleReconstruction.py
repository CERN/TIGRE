"""
A first reconstruction
======================

Reconstruct simulated data with FDK and OS-SART and compare the results.

Source: `Python/demos/d04_SimpleReconstruction.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d04_SimpleReconstruction.py>`__ · `MATLAB/Demos/d04_SimpleReconstruction.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d04_SimpleReconstruction.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d04_SimpleReconstruction.m
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
# Reconstruct image using OS-SART and FDK
# ---------------------------------------

# FDK
imgFDK = algs.fdk(noise_projections, geo, angles)
# OS-SART

niter = 50
imgOSSART = algs.ossart(noise_projections, geo, angles, niter)

# %%
# Show the results
# ----------------

tigre.plotimg(np.concatenate([imgFDK, imgOSSART], axis=1), dim="z")
