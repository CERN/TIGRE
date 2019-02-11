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
geo=defaultGeometry('nVoxel',[128;128;128]);                     

%% Load data and generate projections 
% see previous demo for explanation
angles=linspace(0,2*pi,200);
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
noise_projections=addCTnoise(projections);

%% Usage of FDK

% the FDK algorithm has been taken and modified from 
% 3D Cone beam CT (CBCT) projection backprojection FDK, iterative reconstruction Matlab examples
% https://www.mathworks.com/matlabcentral/fileexchange/35548-3d-cone-beam-ct--cbct--projection-backprojection-fdk--iterative-reconstruction-matlab-examples

% The algorithm takes, as eny of them, 3 mandatory inputs:
% PROJECTIONS: Projection data
% GEOMETRY   : Geometry describing the system
% ANGLES     : Propjection angles
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
imgFDK1=FDK(noise_projections,geo,angles,'filter','hann');
imgFDK2=FDK(noise_projections,geo,angles,'filter','ram-lak');


% They look quite the same
plotImg([imgFDK1 imgFDK2],'Dim','Z');

% but it can be seen that one has bigger errors in the whole image, while
% hte other just in the boundaries
plotImg([abs(head-imgFDK1) abs(head-imgFDK2)],'Dim','Z');
