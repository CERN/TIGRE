%% DEMO 20:  Algorithms 05. FISTA with TV proximal
%
%
%  This demo presents the FISTA algorithm in TIGRE. using a total variation proximal.
%  Total variation algorithms try to minimize the variation (gradient) of the
%  image, assuming its piecewise smooth, as most things in nature are (i.e.
%  human body). FISTA attempts to perform this with much faster convergence
%  than sandard gradient algorithms such as SIRT.
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

%% FISTA is a quadratically converging algorithm.

% 'hyper': This parameter should approximate the largest 
%          eigenvalue in the A matrix in the equation Ax-b and Atb. 
%          Empirical tests show that for, the headphantom object:
%               
%               geo.nVoxel = [64,64,64]'    ,      hyper (approx=) 2.e8
%               geo.nVoxel = [512,512,512]' ,      hyper (approx=) 2.e4
%          Default: 2.e8
% 'tviter':  Number of iterations of Im3ddenoise to use. Default: 20
% 'lambda':  Multiplier for the tvlambda used, which is proportional to 
%            L (hyper). Default: 0.1
% 'verbose': get feedback or not. Default: 1
%--------------------------------------------------------------------------
imgFISTA_default = FISTA(noise_projections,geo,angles,100);

%% Adding a different hyper parameter
imgFISTA_hyper = FISTA(noise_projections,geo,angles,100,'hyper' ,2.e6);

%% Recon more significant tv parameters
imgFISTA_hightv = FISTA(noise_projections,geo,angles,100,'hyper',2.e6,'tviter',100,'lambda',20);

%% Faster converging FISTA
imgFISTA_mod = FISTA_mod(noise_projections,geo,angles,100,'hyper',2.e6,'tviter',100,'lambda',20);

%% Plot results
plotImg([imgFISTA_default,imgFISTA_hyper,imgFISTA_hightv,imgFISTA_mod])