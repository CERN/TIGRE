%% DEMO 14:  Playing with offsets
%
%
% In this demo we show how to change offsets to either the iamge or the
% detector, and the flexibility of it.
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

%% Define Geometry
geo=defaultGeometry('nVoxel',[256;256;256]);                     

% Offsets
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lets try simnple offset: The detector gets completelly displaced
geo.offOrigin = [0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[180; 0];                     % Offset of Detector            (mm)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary 
geo.accuracy=0.5;                           % Accuracy of FWD proj          (vx/sample)

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,360);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');

plotImg(projections, 'Dim',3)

%% lets see it
imgFDK = FDK_half(projections,geo,angles);

plotImg(imgFDK,'Dim',3)
close all
subplot(2,2,1)
imshow(imgFDK(:,:,64), []), title('Slice = 64')
subplot(2,2,2)
imshow(imgFDK(:,:,96), []), title('Slice = 96')
subplot(2,2,3)
imshow(imgFDK(:,:,128), []), title('Slice = 128')
subplot(2,2,4)
imshow(imgFDK(:,:,160), []), title('Slice = 160')
