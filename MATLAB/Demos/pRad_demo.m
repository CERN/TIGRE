%% Creating one CUDA accelerated test radiography

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
% Coded by:           Stefanie Kaser 
%--------------------------------------------------------------------------
%% pCT geometry definition

geo.dDetector = [0.25; 0.25]; % Pixel size of optimized pRad in mm.
geo.DSID = 245; % Distance between source and upstream detector.
geo.DSO = 300; % Distance between source and origin (center of rotation).
geo.DSD = 355; % Distance between source and downstream detector.
geo.hull = [0; 0; 0; 0]; % We won't use a convex hull here (all entries are 
% set to 0).
geo.sDetector = [20; 20]; % Defines the size (in mm) of the optimized pRad.
geo.mode = 'parallel'; % or 'cone' in case of cone beam geometry

%% pCT initial beam energy 
eIn = single(100);

%% Load pCT data set. 
% The data set corresponds to one radiographic image (pRad) of an Aluminum 
% stair profile as used in https://arxiv.org/abs/2106.12890. However, our 
% test data set is based on Monte Carlo simulations (GATE: 
% doi.org/10.1016/S0920-5632(03)90969-8). The data set contains the 
% protons' upstream and downstream positions and directions as well as the 
% single proton WEPLs (all in mm).
data = pCTdata();

%% Binning the data into an optimized pRad
% Finally, the single proton data are binned into an optimized pRad.
proj = pCTCubicSpline_mex(data.posIn, data.posOut, data.dirIn, ...
    data.dirOut, data.Wepl, eIn, geo);
% We are only creating one test pRad here. Creating pRads at multiple 
% rotation angles would allow to use the collection as input for TIGRE's 
% reconstruction algorithms (see d04_SimpleReconstruction.m).

%% Plot the result
imshow(proj, [0, 3], 'InitialMagnification', 500);