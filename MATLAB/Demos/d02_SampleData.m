%% DEMO 02: Sample data in TIGRE
%
%
% TIGRE has some sample data so you can use without the need of having
% your own data. This code sample shows how to load the data.
%
% Sample data is stored in "data" folder, inside TIGRE/Common.
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

geo=defaultGeometry();
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
% head phantom
%
%
head=headPhantom(geo.nVoxel); %default is 128^3
% show it
plotImg(head,'Dim','Z');
citeme('headPhantom')


