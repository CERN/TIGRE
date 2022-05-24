%% DEMO 17: Detector Rotation
%
%
%
%  Some systems have a slight rotation of the detector due to mechanicall
%  inacuracies. 
%
%  According to the article "A geometric calibration method for cone beam
%  CT systems" (DOI: 10.1118/1.2198187), only Roll needs to be corrected
%  for in the algorithmic part, as the other 2 can be easily ignored if
%  "sufficiently small". In TIGRE we decided to leave that to the users
%  discretion and implemented the 3 possible detector rotation, per angle
%  if needed. 
%
%  This demo shows how to use it. 
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

% a bit smaller than usual because the demo includes a very big detector
% angle for showcase
geo.sVoxel=[256;256;256]/2;                 % total size of the image       (mm)


geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)


% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)

%% Geometry has an additional field for the rotation.
%
%
angles=linspace(0,2*pi,100);

% lets define 3 angles, with big variability in each direction

roll=angles;
pitch=0.7*linspace(0,1,size(angles,2));
yaw=0.7*linspace(0,1,size(angles,2));

% Fill rotDetector. It can also be a 3x1, if its not angle dependent. e.g.
% geo.rotDetector=[pi/5; 0; pi/10];

geo.rotDetector=[roll;pitch;yaw];



%% Load data and generate projections 
% define angles
% Load thorax phatom data
head=headPhantom(geo.nVoxel);
% generate projections
projections=Ax(head,geo,angles,'interpolated');

%% Lets plot the projections with a rotated detector
plotProj(projections,angles);

%% lets reconstruct with and without detector rotation
% 
imgRotDet=OS_SART(projections,geo,angles,50);

% No rotation
geo.rotDetector=[0;0;0];
projections2=Ax(head,geo,angles,'interpolated');


imgnoRot=OS_SART(projections2,geo,angles,50);

% %% Plot to show that indeed the reconstruction is right
% 
% 
% 
plotImg([imgRotDet imgnoRot] ,'dim',3)
