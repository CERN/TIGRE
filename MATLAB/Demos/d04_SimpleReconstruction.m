%% Demo 4: Simple Image reconstruction
%
%
% This demo will show how a simple image reconstruction can be performed,
% by using OS-SART and FDK
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------
%% Initialize

clear;
close all;
%% Geometry
geo=defaultGeometry('nVoxel',[128;128;128]);                     

%% Load data and generate projections 
% define angles
angles=linspace(0,2*pi,100);
% Load thorax phatom data
head=headPhantom(geo.nVoxel);
% generate projections
projections=Ax(head,geo,angles,'interpolated');
% add noise
noise_projections=addCTnoise(projections);

%% Reconstruct image using OS-SART and FDK

% FDK
imgFDK=FDK(noise_projections,geo,angles);
% OS-SART

niter=50;
imgOSSART=OS_SART(noise_projections,geo,angles,niter);

% Show the results
plotImg([imgFDK,imgOSSART],'Dim','Z');
