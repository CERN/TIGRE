%% DEMO 9:  Algorithms 06.Anatomical Priors
%
%
%  This demo presents the TIGRE capability of using image priors for
%  CT reconstruction. The most known algorithm for this is PICCS
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
% Coded by:           Ander Biguri, Reuben Lindroos
%--------------------------------------------------------------------------
%% Initialize
clear;
close all;

%% Define Geometry
geo=defaultGeometry('nVoxel',[128;128;128]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi-2*pi/30,30);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);

%% Lets create a OS-SART test for comparison
[imgOSSART,errL2OSSART]=OS_SART(noise_projections,geo,angles,60);

%% ASD_POCS (TV without prior)
epsilon=errL2OSSART(end);
alpha=0.002;
ng=25;

imgASDPOCS=ASD_POCS(noise_projections,geo,angles,50,'TViter',ng,'maxL2err',epsilon,'alpha',alpha);
%% PICCS (ASD_POCS, but with an image prior)
imgPICCS=PICCS(noise_projections,geo,angles,50,head,'TViter',ng,'maxL2err',epsilon,'alpha',alpha);

%% Plot results
errOSSART=abs(imgOSSART-head);
errASDPOCS=abs(imgASDPOCS-head);
errPICCS=abs(imgPICCS-head);

plotImg( [imgOSSART imgASDPOCS imgPICCS],'Dim','z','slice',34);
plotImg( [errOSSART errASDPOCS errPICCS],'Dim','z','slice',34,'clims',[0,0.1]);