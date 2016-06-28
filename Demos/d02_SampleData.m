%% DEMO 02: Sample data in TIGRE
%
%
% TIGRE has some sample data so you can use without the need of having
% your own data. This code sample shows how to load the data.
%
% Sample data is stored in "test_data" folder, inside TIGRE.
%
% If you want to contribute your own phantom/real data, please do. Send an
% email to tigre.toolbox@gmail.com. If you want us to just add a link to
% your own page with data, we could also do that.
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
geo.nVoxel=[256;256;256];                   % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)


% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)

%% Sample data:
%--------------------------------------------------------------------------
%3D Shepp-Logan
%
% Type of shepp logan phantom. The shape will be the same, but Hounsfield
% values will be different
% Options:
%--------------------------------------------------------------------------
% License problems, we are looking at this now.
% shepp_type='yu-ye-wang'; 
% shepp_type='Shepp-Logan'; 
% shepp_type='Modified Shepp-Logan';  % (default).
% 
% shepp=sheppLogan3D(geo.nVoxel,shepp_type); % Default are 128^3 and Modified shepp-logan
%
% Get it here: http://uk.mathworks.com/matlabcentral/fileexchange/50974-3d-shepp-logan-phantom
%
% show it
% plotImg(shepp,'Dim','Z');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Thorax phantom
%
%
thorax=thoraxPhantom(geo.nVoxel); %default is 128^3
% show it
plotImg(thorax,'Dim','Z');



