#%% DEMO 02: Sample data in TIGRE
#
#
# TIGRE has some sample data so you can use without the need of having
# your own data. This code sample shows how to load the data.
#
# Sample data is stored in "data" folder, inside TIGRE/Common.
#
# If you want to contribute your own phantom/real data, please do. Send an
# email to tigre.toolbox@gmail.com. If you want us to just add a link to
# your own page with data, we could also do that.
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
#
import tigre
from tigre.utilities import sample_loader
from tigre.utilities import sl3d

#%% Define geometry
geo = tigre.geometry_default(high_resolution=False)


#%%
# --------------------------------------------------------------------------
# 3D Shepp-Logan
#
# Type of shepp logan phantom. The shape will be the same, but Hounsfield
# values will be different
# Options:
#  phantom_type: str
#      One of {"kak-slaney", "yu-ye-wang", "toft-schabel"}, optional
#      Default is "yu-ye-wang"
#      The type of phantom.
#  size_out : int or list whose length is three, optional
#      Default is [128, 128, 128]
#      The number of voxels of phantom. [nVoxelZ, nVoxelY, nVoxelX]
# --------------------------------------------------------------------------
# phantom_type="kak-slaney"    # Air is 0. Water is 1. Proportional to Hounsfield value.
phantom_type = (
    "yu-ye-wang"  # Default of Python TIGRE Shepp-Logan phantom. Improved visual perception
)
# phantom_type="toft-schabel"  # Default of MATLAB TIGRE Shepp-Logan phantom.
shepp = sl3d.shepp_logan_3d(
    geo.nVoxel, phantom_type=phantom_type
)  # Default are 128^3 and "yu-ye-wang"
tigre.plotImg(shepp, dim="z")

#%% load head phantom
head = sample_loader.load_head_phantom(geo.nVoxel)
# check the shape
print(head.shape)
# show it
tigre.plotImg(head, dim="z")
