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
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------
%% Initialize
clear;
close all;

%% Define Geometry
% 
% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[512; 512];					% number of pixels              (px)
geo.dDetector=[0.8; 0.8]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[128;128;128];                   % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,100);
thorax=thoraxPhantom(geo.nVoxel);
projections=Ax(thorax,geo,angles,'interpolated');
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
plotImg(abs([thorax-imgFDK thorax-imgdenoised thorax-imcroped]),'Dim','Z')


