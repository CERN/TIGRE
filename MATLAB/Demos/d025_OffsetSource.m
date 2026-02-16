%% DEMO 03: Generate sample data and add realistic CT noise to it.
%
% This demo will show how to generate sample data for image reconstruction
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
%% Geometry
geo=defaultGeometry();
head=headPhantom(geo.nVoxel);
angles = linspace(0,2*pi,500);
P0=Ax(head,geo,angles,'interpolated');

img_ref = FDK(P0,geo,angles);
imtool(img_ref(:,:,128),[])

offset = 60;
geo1 = geo;
geo1.offSource = [offset;0];
geo1.offDetector = [0;0];
geo1.offOrigin = [0;0;0];
P1=Ax(head,geo1,angles,'interpolated');

img = FDK(P1,geo1,angles);
% show the difference caused by the offset of source
imtool([P0(:,:,100)-P1(:,:,100)],[])
% show the similar reconstruction images
imtool([img_ref(:,:,128),img(:,:,128)],[])
