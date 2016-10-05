%% Demo 13: Helical Geometry tests
%
%
% This demo shows an example of TIGRE wirking on Helical scan geometries
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
%%

clear;
clc;

geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[512; 512];					% number of pixels              (px)
geo.dDetector=[0.8; 0.8]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[128;128;512];                   % number of voxels              (vx)
geo.sVoxel=[256;256;1024];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
geo.COR=0;
geo.mode='cone';

% Auxiliary
geo.accuracy=0.5;       


%% Create data+ angles
angles=linspace(00.0,2*pi,100);
angles=[angles angles angles];
% Load thorax phatom data
thorax=thoraxPhantom(geo.nVoxel); % yes, not the best example data, but It will do.

% Thsi makes it helical
geo.offOrigin(3,:)=linspace(-1024/2+128,1024/2-128,length(angles)).';  % about 256^3 images fit int he detector with this size. 
geo.offOrigin(1,:)=0;
geo.offOrigin(2,:)=0;

% project data
data=Ax(thorax,geo,angles);

% Uncomment if you want to see the data
% plotProj(data,angles);
%% Reconstruct Helical

OSSARTimg=OS_SART(data,geo,angles,30);
% SARTimg=SART(data,geo,angles,30); % takes time
CGLSimg=CGLS(data,geo,angles,20);
%% Plot results

% CGLS and SIRT
plotImg([thorax, OSSARTimg ,CGLSimg],'Dim',3,'Step',3);


