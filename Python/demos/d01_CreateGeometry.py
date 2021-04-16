#%% DEMO 01: Describing your geometry
#
#  In TIGRE the geometry is stored in an structure. To see documentation
#  about geometry, run:
#
#     doc('TIGRE/Geometry')
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
#%% Import statements
import tigre
import numpy as np

#%%  Geometry definition
#
#                  Detector plane, behind
#              |-----------------------------|
#              |                             |
#              |                             |
#              |                             |
#  Centered    |                             |
#    at O      A V    +--------+             |
#              |     /        /|             |
#     A Z      |    /        / |*D           |
#     |        |   +--------+  |             |
#     |        |   |        |  |             |
#     |        |   |     *O |  +             |
#     *--->y   |   |        | /              |
#    /         |   |        |/               |
#   V X        |   +--------+        U       |
#              .--------------------->-------|
#
#            *S
#
#
# An image of the geometry can be found at DOI: 10.1088/2057-1976/2/5/055010
# or simply at: https://i.imgur.com/mRweux3.png
#%% Geometry class:
#           -nVoxel:        3x1 array of number of voxels in the image
#           -sVoxel:        3x1 array with the total size in mm of the image
#           -dVoxel:        3x1 array with the size of each of the voxels in mm
#           -nDetector:     2x1 array of number of voxels in the detector plane
#           -sDetector:     2x1 array with the total size in mm of the detector
#           -dDetector:     2x1 array with the size of each of the pixels in the detector in mm
#           -DSD:           1x1 or 1xN array. Distance Source Detector, in mm
#           -DSO:           1x1 or 1xN array. Distance Source Origin.
#           -offOrigin:     3x1 or 3xN array with the offset in mm of the centre of the image from the origin.
#           -offDetector:   2x1 or 2xN array with the offset in mm of the centre of the detector from the x axis
#           -rotDetector:   3x1 or 3xN array with the rotation in roll-pitch-yaw of the detector
#%% Example
#
# Lets create a geomtry object
geo = tigre.geometry()
# VARIABLE                                   DESCRIPTION                    UNITS
# -------------------------------------------------------------------------------------
# Distances
geo.DSD = 1536  # Distance Source Detector      (mm)
geo.DSO = 1000  # Distance Source Origin        (mm)
# Detector parameters
geo.nDetector = np.array([512, 512])  # number of pixels              (px)
geo.dDetector = np.array([0.8, 0.8])  # size of each pixel            (mm)
geo.sDetector = geo.nDetector * geo.dDetector  # total size of the detector    (mm)
# Image parameters
geo.nVoxel = np.array([256, 256, 256])  # number of voxels              (vx)
geo.sVoxel = np.array([256, 256, 256])  # total size of the image       (mm)
geo.dVoxel = geo.sVoxel / geo.nVoxel  # size of each voxel            (mm)
# Offsets
geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin   (mm)
geo.offDetector = np.array([0, 0])  # Offset of Detector            (mm)
# These two can be also defined
# per angle

# Auxiliary
geo.accuracy = 0.5  # Variable to define accuracy of
# 'interpolated' projection
# It defines the amoutn of
# samples per voxel.
# Recommended <=0.5             (vx/sample)

# Optional Parameters
# There is no need to define these unless you actually need them in your
# reconstruction


geo.COR = 0  # y direction displacement for
# centre of rotation
# correction                   (mm)
# This can also be defined per
# angle

geo.rotDetector = np.array([0, 0, 0])  # Rotation of the detector, by
# X,Y and Z axis respectively. (rad)
# This can also be defined per
# angle

geo.mode = "cone"  # Or 'parallel'. Geometry type.

#%% Print the Geometry
print(geo)
#%% Alternatively,
# # if you are just experimenting with TIGRE and you dont care too much about the geometry, you can generate default geometries as:
geo = tigre.geometry_default()  # Default cone beam geometry
geo = tigre.geometry(
    mode="cone", default=True
)  # default=True calculates all parameters for the geometry for you, so you can do:
geo = tigre.geometry(
    mode="cone", default=True, nVoxel=np.array([256, 256, 256])
)  # This will calculate a reasonable geometry for this number of voxels
geo = tigre.geometry(
    mode="cone", default=True, high_resolution=True
)  # high_resolution will result on an image 512^3, while false will result on a 128^3 image. nVoxel overrides these values.
geo = tigre.geometry(
    mode="parallel", nVoxel=np.array([512, 512, 512])
)  # Parallel beam geometry does not require anything other than the image size.

#%% Plot your geometry
geo = tigre.geometry_default()  # Default cone beam geometry
tigre.plot_geometry(geo, angle=-np.pi / 6)
