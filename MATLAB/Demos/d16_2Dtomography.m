%% DEMO 16: TIGRE 2D tomography 
%
%
%  In demo we show how to do 2D tomography in TIGRE. It is exactly the same
%  as 3D.
%
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
%% PARALLEL BEAM 2D

% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters

%
% 1- Make the detector Nx1 
% 2- Make the image NxMx1 


% Image parameters
geo.nVoxel=[256;256;1];                   % number of voxels              (vx)
geo.sVoxel=[256;256;1];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)

geo.nDetector=[512;  1];					% number of pixels              (px)
geo.dDetector=[0.8; geo.dVoxel(3)]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)

% MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!

% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0]; 


% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)
geo.mode='parallel';
%% Define angles of projection and load phatom image

angles=linspace(0,2*pi,100);
phatom=single(phantom('Modified Shepp-Logan',geo.nVoxel(1)));
projections=Ax(phatom,geo,angles);
%% recosntruct

imgOSSART=OS_SART(projections,geo,angles,40);
imgASDPOCS=ASD_POCS(projections,geo,angles,40);

%% Plot
plotImg([imgOSSART;imgASDPOCS ],'Dim',3,'Slice',1); %top is ASD_POCS, bottom OSSART


%%
%%
%% And now Fan Beam
%% 
% The same thing!

%% FAN BEAM 2D

% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters

% Image parameters
geo.nVoxel=[256;256;1];                   % number of voxels              (vx)
geo.sVoxel=[256;256;1];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)

geo.nDetector=[512;  1];					% number of pixels              (px)
geo.dDetector=[0.8; 1]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)

% Offsets
geo.offOrigin =[0;0;0];                        % Offset of image from origin   (mm)              ; 
geo.offDetector=[0; 0]; 

% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)
geo.mode='cone';
%% Define angles of projection and load phatom image

angles=linspace(0,2*pi,100);
phatom=single(phantom('Modified Shepp-Logan',geo.nVoxel(1)));
% phatom=cat(3,phatom,phatom);
projections=Ax(phatom,geo,angles,'interpolated');
%% recosntruct

imgOSSART=OS_SART(projections,geo,angles,40);
% in 2D, the TV norm minimization happens in CPU, this causes the
% algorithms to be slower
imgASDPOCS=ASD_POCS(projections,geo,angles,40);

%% Plot
plotImg([imgOSSART;imgASDPOCS ],'Dim',3,'Slice',1); %top is ASD_POCS, nottom OSSART


