%% DEMO 03: Generate sample data and add realistic CT noise to it.
%
% This demo will show how to generate sample data for image reconstruction
%
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
geo=defaultGeometry();
%% Define angles of projection and load phatom image

% define projection angles (in radians)
angles=linspace(0,2*pi,100);
% load phatnom image
head=headPhantom(geo.nVoxel);

% Simulate forward projection.
% Strongly suggested to use 'iterpolated' option for more accurate
% projections. reduce geo.accuracy for better results
projections=Ax(head,geo,angles,'interpolated');

% Add realistic noise. Adds photon scattering noise ('Poisson') and
% electronic noise of the detector ('Gaussian').
%
% 'Poisson' is related to the maximum photon count in the detector. 1e5 is
% a standard clinical nuber, reduce it for more noise
% 'Gaussian' is related to possible electronic noise in the detector. mean
% of 0 and std of 10 is common in clinical scenario. Increase std for more
% noise.
noise_projections=addCTnoise(projections,'Poisson',1e5,'Gaussian',[0 10]);

% Plot Projections
plotProj(projections,angles)
% plot noise
plotProj(projections-noise_projections,angles)


