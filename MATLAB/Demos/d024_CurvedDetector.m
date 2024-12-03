%% Demo 4: Simple Image reconstruction
%
%
% This demo will show how a simple image reconstruction can be performed,
% by using OS-SART and FDK
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
geo=defaultGeometry('nVoxel',[128;128;128]);                     

% If you have a curved detector, currently the only way to process it is to "flatten" the detector. 

proj= ...
% You can change the oversample value to create denser (and thus more accurate) flat projections, in exchange for computational time at recon
oversample = 1; % this is optional 
proj = flatten_detector(proj,geo,oversample);


% Now just use as if proj came from flat panel detectors. 