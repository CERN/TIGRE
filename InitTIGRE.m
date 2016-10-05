% This Function initializes the toolbox.
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

% Add tolbox folders
addpath('./Algorithms');
addpath('./Utilities');
addpath('./Utilities/Quality_measures');
addpath(genpath('./Test_data'));

% different arch versions
if ispc
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./Mex_files/win64');
    else
        addpath('./Mex_files/win32');
    end
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./Mex_files/mac64');
    else
        addpath('./Mex_files/mac32');
    end
else
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./Mex_files/linux64');
    else
        addpath('./Mex_files/linux32');
    end
end
    
addpath('./Demos');

% Perceptually uniform colormaps
addpath('./Colormaps');
% Add third party tools from FEX
addpath('./Third_party_tools/arrow3d'); % 3D shepp-Logan
addpath('./Third_party_tools/sec2hours');
addpath('./Third_party_tools/readMHD');

if ispc
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./Mex_files/win64');
    else
        addpath('./Mex_files/win32');
    end
elseif ismac
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./Mex_files/mac64');
    else
        addpath('./Mex_files/mac32');
    end
else
    if ~isempty(strfind(computer('arch'),'64'))
        addpath('./Mex_files/linux64');
    else
        addpath('./Mex_files/linux32');
    end
end

if ispc
    [user, sys]=memory;
    
    if sys.PhysicalMemory.Total<9000000000 % 8Gb
        warning('Your Computer has 8Gb or less of RAM memory. Using image sizes of higher than 512^3 is not recomended (most likely not possible)')
    end
    
    if sys.PhysicalMemory.Total<2500000000 % 2Gb
        warning('Your Computer has 2Gb or less of RAM memory. Using image sizes of higher than 256^3 is not recomended (most likely not possible)')
    end
else
    warning('TIGRE needs a big amount of memory, be careful when running big images.')
end

clear all;
