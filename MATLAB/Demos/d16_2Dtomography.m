%% DEMO 16: Tricking TIGRE to do 2D tomography 
%
%
%  In this more advanced demo, a way of tricking TIGRE to do 2D tomography
%  is shown. 
%
%  TIGRE by default does not support 2D tomography, as the GPU acceleration
%  is not that relevant when the problem is smaller (as in 2D). However,
%  one may be interested in testing the algorithms included in TIGRE for 2D
%  cases, thus a way of doing it is shown. 
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

% HERE COMES THE TRICK.
%
% 1- Make the detector Nx2 (instead of 1)
% 2- Make the image NxMx2  (isntead of 1)


% Image parameters
geo.nVoxel=[256;256;2];                   % number of voxels              (vx)
geo.sVoxel=[256;256;2];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)

geo.nDetector=[512;  2];					% number of pixels              (px)
geo.dDetector=[0.8; geo.dVoxel(3)]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)

% MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!

% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)


% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)
geo.mode='parallel';
%% Define angles of projection and load phatom image

angles=linspace(0,2*pi,100);
phatom=single(phantom('Modified Shepp-Logan',geo.nVoxel(1)));
phatom=cat(3,phatom,phatom);
projections=Ax(phatom,geo,angles,'interpolated');
%% recosntruct

imgOSSART=OS_SART(projections,geo,angles,40);
imgASDPOCS=ASD_POCS(projections,geo,angles,40);

%% Plot
plotImg([imgOSSART;imgASDPOCS ],'Dim',3,'Slice',1); %top is ASD_POCS, nottom OSSART


%%
%%
%% And now Fan Beam
%% 
% This will be considerably trickier, but does work anyways. As we want to
% use the cone beam geometry to simulate fan beam, we need to make sure
% that the second row of the detector will never be used in the
% backprojection. This means:
% 1- We need to offset the detector in Z so it goes away
% 2- The detector pixel size has to be VERY big in V, big enough to
% completelly miss the image in the backprojection

%% FAN BEAM 2D

% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters

% Image parameters
geo.nVoxel=[256;256;2];                   % number of voxels              (vx)
geo.sVoxel=[256;256;2];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)

% This assumes 0 offsets in the image
vsize=(geo.dVoxel(3)*geo.DSD)/(geo.DSO-sqrt((geo.sVoxel(1)/2).^2+(geo.sVoxel(2)/2).^2));
vsize=2*vsize;

geo.nDetector=[512;  2];					% number of pixels              (px)
geo.dDetector=[0.8; vsize]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)

% MAKE SURE THAT THE DETECTOR PIXELS SIZE IN V IS THE SAME AS THE IMAGE!

% Offsets
geo.offOrigin =[0;0;0];                        % Offset of image from origin   (mm)              
geo.offDetector=[0; -vsize/2];                  % Offset of Detector            (mm)


% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)
geo.mode='cone';
%% Define angles of projection and load phatom image

angles=linspace(0,2*pi,100);
phatom=single(phantom('Modified Shepp-Logan',geo.nVoxel(1)));
phatom=cat(3,phatom,phatom);
projections=Ax(phatom,geo,angles,'interpolated');
%% recosntruct

imgOSSART=OS_SART(projections,geo,angles,40);
imgASDPOCS=ASD_POCS(projections,geo,angles,40);

%% Plot
plotImg([imgOSSART;imgASDPOCS ],'Dim',3,'Slice',1); %top is ASD_POCS, nottom OSSART


