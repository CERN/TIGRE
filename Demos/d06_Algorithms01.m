%% DEMO 6:  Algorithms 01.
%
%
% In this demo the usage of the FDK is explained
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

%% Usage of FDK

% the FDK algorithm has been taken and modified from 
% 3D Cone beam CT (CBCT) projection backprojection FDK, iterative reconstruction Matlab examples
% https://www.mathworks.com/matlabcentral/fileexchange/35548-3d-cone-beam-ct--cbct--projection-backprojection-fdk--iterative-reconstruction-matlab-examples

% The algorithm takes, as eny of them, 3 mandatory inputs:
% PROJECTIONS: Projection data
% GEOMETRY   : Geometry describing the system
% AMGLES     : Propjection angles
% And has a single optional argument:
% FILTER: filter type applied to the projections. Possible options are
%        'ram-lal' (default)
%        'shepp-logan'
%        'cosine'
%        'hamming'
%        'hann'
% The choice of filter will modify the noise and sopme discreatization
% errors, depending on which is chosen.
%
imgFDK1=FDK(noise_projections,geo,angles,'hann');
imgFDK2=FDK(noise_projections,geo,angles,'ram-lak');


% They look quite the same
plotImg([imgFDK1 imgFDK2],'Dim','Z');

% but it can be seen that one has bigger errors in the whole image, while
% hte other just in the boundaries
plotImg([abs(thorax-imgFDK1) abs(thorax-imgFDK2)],'Dim','Z');
