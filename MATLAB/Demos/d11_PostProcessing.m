%% Demo 11: Postprocessing 
%
% This demo demonstrates the available postprocessing tools in TIGRE.
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

%% Define Geometry
geo=defaultGeometry('nVoxel',[128;128;128]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,100);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);
%% Lets just use FDK

imgFDK=FDK(noise_projections,geo,angles);

%% Postprocessing
%
%  Currently the postprocessing steps in TIGRE are limited, but there are 2
%  fucntions you can use
%
% im3Ddenoise : Denoises a 3D image, using Total Variation denoising.
%
%  Argumetns are the hyperparameter and number of iterations (same as in
%  SART-TV)
imgdenoised=im3DDenoise(imgFDK,'TV',100,15);

% cropCBCT: Crops all the sections that lie outide the area that is covered
% by the cone
%
%
imcroped=cropCBCT(imgFDK,geo);

%% plot results
% denoised image is clearer
plotImg([imgFDK imgdenoised imcroped],'Dim','Z')
% however, teh denoising has no knoledge of the original data (projections)
% this it doesnt reduce the error. The error increases, specially in small
% areas
plotImg(abs([head-imgFDK head-imgdenoised head-imcroped]),'Dim','Z')


