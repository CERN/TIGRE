##% Demo 5: How to use plotting functions
#
# This demo will demostrate the options for plotting projection and images
# on TIGRE. The functions have been in previous demos, but in here an
# exaustive explanation and usage of them is given.
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
#%%Initialize
import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs

#%% Geometry
geo=tigre.geometry_default(high_resolution=False) 

#%% Load data and generate projections 
# define angles
angles=np.linspace(0,2*np.pi,100)
# Load thorax phatom data
head = sample_loader.load_head_phantom(geo.nVoxel)
# generate projections
projections=tigre.Ax(head,geo,angles)
# add noise
noise_projections=CTnoise.add(projections,Poisson=1e5,Gaussian=np.array([0, 10]))

#%% Reconstruct image using OS-SART and FDK

# FDK
imgFDK=algs.fdk(noise_projections,geo,angles)
niter=50
imgOSSART=algs.ossart(noise_projections,geo,angles,50)

#%% %% Lets use PlotProj

# See issue #265


#%% PlotImg

# plotImg plots the image slice by slice. 
#
# List of optional parameters: 


# 'Dim': specifies the dimension for plotting. 
#    Dim can be 'X','Y','Z'
#
dimension='Z' 

# 'Step': step size of the plotting. Useful when images are big or one just
# wants an overview of the result

step=2
# 'Colormap': Defines the colormap used to plot. Default is 'gray'. 

colormap='plasma'
colormap='magma'
colormap='gray'
colormap='viridis'

# 'Clims': Defines the data limits for the color, usefull when one wants to
# see some specific range. The default computes the 5% and 95% percentiles
# of the data and uses that as limit.

clims=[0, 1]

# 'Savegif': allows to save the plotted figure as an animated gif,
# specified by the given filename.

giffilename='demo5image.gif'

# 'Slice': allows to plot a single slice .Will overwrite the behaviour
# of 'Step'
slice=64

# Lets go for it

tigre.plotimg(imgFDK,dim=dimension,step=step,clims=clims,colormap=colormap,savegif=giffilename)

# Remember: You can always plot more than 1 image together!
tigre.plotimg(np.concatenate([head, imgFDK, imgOSSART],axis=1),dim='z')
# Or even the errors!
plotImg(np.concatenate([np.abs(head-imgFDK), np.abs(head-imgOSSART)]),dim='z')