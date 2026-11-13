"""
Arbitrary axis of rotation
==========================

Trajectories that are not a circle around z, given as Euler angles.

Source: `Python/demos/d18_ArbitraryAxisOfRotation.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d18_ArbitraryAxisOfRotation.py>`__ · `MATLAB/Demos/d18_ArbitraryAxisOfRotation.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d18_ArbitraryAxisOfRotation.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d18_ArbitraryAxisOfRotation.m
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

## Define angles
numProjs = 100

anglesY = np.linspace(0, 2 * np.pi, numProjs, endpoint=False)
anglesZ2 = anglesY
anglesZ1 = np.pi * np.sin(anglesY)
angles = np.vstack([anglesZ1, anglesY, anglesZ2]).T

## Get Image

head = sample_loader.load_head_phantom(geo.nVoxel)

## Project

projections = tigre.Ax(head, geo, angles)

tigre.plotproj(projections, anglesY)  # angle information not right in the title
## Reconstruct:

# Note, FDK will not work.

imgSIRT = algs.sirt(projections, geo, angles, 50)
imgCGLS = algs.cgls(projections, geo, angles, 10)

tigre.plotimg(np.concatenate([head, imgSIRT, imgCGLS], axis=1), dim="z")
