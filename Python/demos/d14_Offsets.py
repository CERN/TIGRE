#%% DEMO 14:  Playing with offsets
#
#
# In this demo we show how to change offsets to either the iamge or the
# detector, and the flexibility of it.
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
## Initialize
#%%Initialize
import tigre
import numpy as np
from tigre.utilities import sample_loader
import tigre.algorithms as algs

#%% Geometry
geo = tigre.geometry_default(high_resolution=False)


# Offsets
## ####################################################################################
# Lets try simnple offset: The detector gets completelly displaced
geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
geo.offDetector = np.array([200, 200])  # Offset of Detector            (mm)
## ####################################################################################
# Auxiliary
geo.accuracy = 0.5
# Accuracy of FWD proj          (vx/sample)

## Load data and generate projections
# see previous demo for explanation
angles = np.linspace(0, 2 * np.pi, 100)

# Load thorax phatom data
head = sample_loader.load_head_phantom(geo.nVoxel)
projections = tigre.Ax(head, geo, angles)


#%% lets see it
tigre.plotproj(projections, angles)

## we will skip reconstruction of this tests because the image is outside the detector


## #####################################################################

#%% Second test: lets test variying offsets:

geo.offDetector = np.vstack(
    [10 * np.sin(angles), 20 * np.cos(angles)]
).T  # Offset of Detector            (mm)
projections2 = tigre.Ax(head, geo, angles)
## lets see it
tigre.plotproj(projections2, angles)
## reconstruction
res = algs.sart(projections2, geo, angles, 10)
tigre.plotimg(res, dim="z")

#%% Third test: lets vary everything

# Lets make the image smaller
geo.nVoxel = np.array([128, 128, 128])  # number of voxels              (vx)
geo.sVoxel = np.array([256, 256, 256]) / 2  # total size of the image       (mm)
geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
head = sample_loader.load_head_phantom(geo.nVoxel)

geo.offDetector = np.vstack(
    [10 * np.sin(angles), 10 * np.cos(angles)]
).T  # Offset of Detector            (mm)
geo.offOrigin = np.vstack(
    [40 * np.sin(angles), np.linspace(-30, 30, 100), 40 * np.cos(angles)]
).T  # Offset of image from origin   (mm)

projections3 = tigre.Ax(head, geo, angles)
## lets see it
tigre.plotproj(projections3, angles)
## reconstruction
res = algs.sart(projections3, geo, angles, 10)
tigre.plotimg(res, dim="z")
