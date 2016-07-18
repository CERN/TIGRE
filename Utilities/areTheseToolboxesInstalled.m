function tf = areTheseToolboxesInstalled(requiredToolboxes)
%ARETHESETOOLBOXESINSTALLED takes a cell array of toolbox names and checks whether they are currently installed
% SYNOPSIS tf = areTheseToolboxesInstalled(requiredToolboxes)
%
% INPUT requiredToolboxes: cell array with toolbox names to test for. Eg. 
%        {'MATLAB','Image Processing Toolbox'}
%
% OUTPUT tf: true or false if the required toolboxes are installed or not
%--------------------------------------------------------------------------
% http://stackoverflow.com/a/2061075/1485872
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

% get all installed toolbox names
v = ver;
% collect the names in a cell array
[installedToolboxes{1:length(v)}] = deal(v.Name);

% check 
tf = all(ismember(requiredToolboxes,installedToolboxes));