"""
Detector rotation
=================

Describe a detector that is rotated or tilted.

Source: `Python/demos/d17_DetectorRotation.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d17_DetectorRotation.py>`__ · `MATLAB/Demos/d17_DetectorRotation.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d17_DetectorRotation.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d17_DetectorRotation.m
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

roll = angles
pitch = 0.7 * np.linspace(0, 1, len(angles))
yaw = 0.7 * np.linspace(0, 1, len(angles))

geo.rotDetector = np.array(np.vstack([roll, pitch, yaw])).T

# %%

# Load thorax phantom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections = tigre.Ax(head, geo, angles)
# add noise
noise_projections = CTnoise.add(projections, Poisson=1e5, Gaussian=np.array([0, 10]))

tigre.plotproj(projections, angles)

# %%
# recon
# -----

imgRotDet = algs.ossart(projections, geo, angles, 50)

#% No rotation
geo.rotDetector = np.array([0, 0, 0])
projections2 = tigre.Ax(head, geo, angles)

imgnoRot = algs.ossart(projections2, geo, angles, 50)

# %%
# Plot to show that indeed the reconstruction is right
# ----------------------------------------------------

tigre.plotimg(np.concatenate([imgRotDet, imgnoRot], axis=1), dim="z")
