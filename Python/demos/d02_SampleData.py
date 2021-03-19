## DEMO 02: Sample data in TIGRE
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
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
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
#--------------------------------------------------------------------------
#
import tigre
from tigre.utilities import sample_loader

#%% Define geometry
geo=tigre.geometry_default() 
#%% load head phantom
head = sample_loader.load_head_phantom(geo.nVoxel)
# show it
tigre.plotImg(head,'Dim','Z')