## Demo 21: Loading scanner data to TIGRE. 
#
# This demo will demostrate the options for loading scanner data into
# TIGRE.
#
#   Supported Manufacturers:
#
#       (not yet, work in progress) Varian     
#
#       Nikon
#
#       Xradia (Zeiss) 
#
#       Bruker Skyscan
#
#       (not yet, work in progress) Philips (in theory, any DICOM, but only tested in Philips Allura) 
#
# Currently we have instructions for generic, Nikon (micro-CT) and Varian
# scanners. 
#
# We are always looking to expand this, if you have code for a scanner that
# is not supported by TIGRE and are allowed to share it, please do, we'd
# like to have as many as possible. 

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
##
import tigre.utilities.io as tigreio
import tigre.algorithms as algs
import tigre
#%% Varian onboard CBCT.

###### TODO

#%% Nikon micro-CT

# Nikon dataset can be loaded with the following code:


datafolder='~/your_data_path/Nikon/Sample_name/'
proj,geo, angles  = tigreio.NikonDataLoader(datafolder)

# as micro-CT datasets are large, optional arguments for loading partial
# amount of data is available:

# load equidistant angles, but only few:
proj, geo, angles  = tigreio.NikonDataLoader(datafolder,sampling='equidistant',num_angles=150)
# load every X angles (10)
proj, geo, angles  = tigreio.NikonDataLoader(datafolder,sampling='step',sampling_step=10)
# load first X angles (1000)
proj, geo, angles  = tigreio.NikonDataLoader(datafolder,sampling='continuous',num_angles=1000)

# You can directly call reconstruction code now:

img=algs.ossart(proj,geo,angles,100)
img=algs.fdk(proj,geo,angles)


#%% Brucker Skyscan

datafolder='~/your_data_path/Bruker/Sample_name/'
proj,geo, angles  = tigreio.BrukerDataLoader(datafolder)

# the same options for sampling and number of angles that exist for Nikon (avobe) exist for Bruker data loaders. 

# Sometimes a folder will contain more than one dataset. 
proj,geo, angles  = tigreio.BrukerDataLoader(datafolder,dataset_num='all') # load all of them
proj,geo, angles  = tigreio.BrukerDataLoader(datafolder,dataset_num=0) # load the first


#%% YXLON 
datafolder='~/your_data_path/YXLON/Sample_name/'
proj,geo, angles  = tigreio.YXLONDataLoader(datafolder)

# the same options for sampling and number of angles that exist for Nikon (avobe) exist for Bruker data loaders. 


#%% DICOM data (only tested on Philips Allura)

######## TODO

#%% Generic

# It is possible that your scanner is not currently supported, or that it
# simply does not have any way of storing the information (Zeiss Xradia
# does not store anything but the projecions)

# if this is the case, the general way of tackling the problem is:

# Step 1: define your geometry and angles
geo=tigre.geometry()
angles=

# Step 2: Load projections:
#store in NP array. Use whichever python lirbary will help you get the data loaded. 


# Step 3: validate

tigre.plotproj(proj,angles)
# you need to make sure that:
#     1) white=metal/bone/high density and black=air.
#         If this is not the case
proj=-np.log(proj/(np.max(proj+1)); # Beer-Lanbert law
#     2) rotation happens left-right instead of top-bottom.
#        If its top bottom permute your data as required

# Step 4: test
imgfdk=algs.fdk(proj,geo,angles)
tigre.plotimg(imgFDK,dim='z')
# If this does not look good, possible things to check:
# - Are the angles in the right direction? maybe they need to be inverted. 
# - Is your geometry properly defined? mistake on geometry will be
#   detrimental to image quality
# - if microCT: halos around the features? COR correction needed. 
