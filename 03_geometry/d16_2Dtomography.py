"""
2D tomography
=============

Fan-beam and parallel-beam 2D problems, as one-slice 3D problems.

Source: `Python/demos/d16_2Dtomography.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d16_2Dtomography.py>`__ · `MATLAB/Demos/d16_2Dtomography.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d16_2Dtomography.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d16_2Dtomography.m
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
# 2D Parallel beam
# ----------------

# VARIABLE                                   DESCRIPTION                    UNITS

# HERE COMES THE TRICK.
#
# 1- Make the detector Nx1 (instead of 1)
# 2- Make the image NxMx1  (instead of 1)

geo = tigre.geometry()
# VARIABLE                                   DESCRIPTION                    UNITS
# -------------------------------------------------------------------------------------
# Distances
geo.DSD = 1536  # Distance Source Detector      (mm)
geo.DSO = 1000  # Distance Source Origin        (mm)
# Image parameters
geo.nVoxel = np.array([1, 256, 256])  # number of voxels              (vx)
geo.sVoxel = np.array([1, 256, 256])  # total size of the image       (mm)
geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
# Detector parameters
geo.nDetector = np.array([1, 512])  # number of pixels              (px)
geo.dDetector = np.array([geo.dVoxel[0], 0.8])  # size of each pixel            (mm)
geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm)
# Offsets
geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
# MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!

geo.mode = "parallel"

# %%
# Define angles of projection and load phantom image
# --------------------------------------------------

angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
head = sample_loader.load_head_phantom(geo.nVoxel)
projections = tigre.Ax(head, geo, angles)
tigre.plotSinogram(projections, 0)

# %%
# reconstruct
# -----------

imgOSSART = algs.ossart(projections, geo, angles, 40)
imgCGLS = algs.cgls(projections, geo, angles, 40)

# %%
# Plot
# ----

tigre.plotimg(np.concatenate([head, imgOSSART, imgCGLS], axis=1), dim="z", slice=0, clims=[0, 0.5])

# %%
# And now Fan Beam
# ----------------

##
# The same thing! Just remember to offset the detector half a pixel, so its
# centered (this is a bug, it will be changed at some point)

## FAN BEAM 2D
geo = tigre.geometry()
# VARIABLE                                   DESCRIPTION                    UNITS
# -------------------------------------------------------------------------------------
# Distances
geo.DSD = 1536  # Distance Source Detector      (mm)
geo.DSO = 1000  # Distance Source Origin        (mm)
# Image parameters
geo.nVoxel = np.array([1, 256, 256])  # number of voxels              (vx)
geo.sVoxel = np.array([1, 256, 256])  # total size of the image       (mm)
geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
# Detector parameters
geo.nDetector = np.array([1, 512])  # number of pixels              (px)
geo.dDetector = np.array([geo.dVoxel[0], 0.8])  # size of each pixel            (mm)
geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm)
# Offsets
geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
# MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!

geo.mode = "cone"

# %%
# Define angles of projection and load phantom image
# --------------------------------------------------

angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
head = sample_loader.load_head_phantom(geo.nVoxel)
projections = tigre.Ax(head, geo, angles)

# %%
# reconstruct
# -----------

imgOSSART = algs.ossart(projections, geo, angles, 40)
imgCGLS = algs.cgls(projections, geo, angles, 40)

# %%
# Plot
# ----

tigre.plotimg(np.concatenate([head, imgOSSART, imgCGLS], axis=1), dim="z", slice=0, clims=[0, 0.5])
