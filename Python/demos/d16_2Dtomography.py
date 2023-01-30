## DEMO 16:  TIGRE 2D tomography
#
#
#  In demo we show how to do 2D tomography in TIGRE. It is exactly the same
#  as 3D.
#
#
#
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

#%% 2D Parallel beam
# VARIABLE                                   DESCRIPTION                    UNITS

# HERE COMES THE TRICK.
#
# 1- Make the detector Nx1 (instead of 1)
# 2- Make the image NxMx1  (isntead of 1)

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

#%% Define angles of projection and load phatom image

angles = np.linspace(0, 2 * np.pi, 100)
head = sample_loader.load_head_phantom(geo.nVoxel)
projections = tigre.Ax(head, geo, angles)

tigre.plotSinogram(projections, 0)
#%% recosntruct

imgOSSART = algs.ossart(projections, geo, angles, 40)
imgCGLS = algs.cgls(projections, geo, angles, 40)

#%% Plot
tigre.plotimg(np.concatenate([head, imgOSSART, imgCGLS], axis=1), dim="z", slice=0, clims=[0, 0.5])


#%% And now Fan Beam
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

#%% Define angles of projection and load phatom image

angles = np.linspace(0, 2 * np.pi, 100)
head = sample_loader.load_head_phantom(geo.nVoxel)
projections = tigre.Ax(head, geo, angles)
#%% recosntruct

imgOSSART = algs.ossart(projections, geo, angles, 40)
imgCGLS = algs.cgls(projections, geo, angles, 40)

#%% Plot
tigre.plotimg(np.concatenate([head, imgOSSART, imgCGLS], axis=1), dim="z", slice=0, clims=[0, 0.5])
