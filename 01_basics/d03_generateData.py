"""
Simulate projections and noise
==============================

Forward-project an image and add realistic Poisson and Gaussian noise.

Source: `Python/demos/d03_generateData.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d03_generateData.py>`__ · `MATLAB/Demos/d03_generateData.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d03_generateData.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d03_generateData.m
      :language: matlab
"""

# %%

import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise

geo = tigre.geometry_default(high_resolution=False)

# %%
# Define angles of projection and load phantom image
# --------------------------------------------------

# define projection angles (in radians)
angles = np.linspace(0, 2 * np.pi, 50, endpoint=False)
# load phantom image
head = sample_loader.load_head_phantom(geo.nVoxel)

# Simulate forward projection.
# To match with mathematical notation, the projection operation is called Ax
projections = tigre.Ax(head, geo, angles)

# Add realistic noise. Adds photon scattering noise ('Poisson') and
# electronic noise of the detector ('Gaussian').
#
# 'Poisson' is related to the maximum photon count in the detector. 1e5 is
# a standard clinical number, reduce it for more noise
# 'Gaussian' is related to possible electronic noise in the detector. mean
# of 0 and std of 10 is common in clinical scenario. Increase std for more
# noise.
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

# %%
# Plot Projections
# ----------------

tigre.plotproj(projections)
# plot noise
tigre.plotproj(projections - noise_projections)
