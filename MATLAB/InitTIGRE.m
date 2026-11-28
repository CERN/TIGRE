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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------

% Add tolbox folders
addpath('./Algorithms');
addpath('./Utilities');
addpath('./Utilities/Quality_measures');
addpath('./Utilities/IO/VarianCBCT');
addpath('./Utilities/IO/Nikon');
addpath('./Utilities/IO/Dicom');
addpath('./Utilities/IO/Bruker');
addpath('./Utilities/IO/YXLON');
addpath('./Utilities/IO/Diondo');
addpath('./Utilities/IO/DXChange');
addpath('./Utilities/GPU');
addpath('./Utilities/DE');

addpath(genpath('./Test_data'));

    
addpath('./Demos');

% Perceptually uniform colormaps
addpath('./Colormaps');
% Add third party tools from FEX
addpath('./Third_party_tools/arrow3d'); % 3D shepp-Logan
addpath('./Third_party_tools/sec2hours');
addpath('./Third_party_tools/readMHD');
addpath('./Third_party_tools/Inpaint_nans/');

% Compiled mex files, added last so they come first on the path
if ispc
    archdir='win';
elseif ismac
    archdir='mac';
else
    archdir='linux';
end
if ~isempty(strfind(computer('arch'),'64'))
    archdir=[archdir '64'];
else
    archdir=[archdir '32'];
end
addpath(['./Mex_files/' archdir]);
% Only there once CompilePCT has been run (proton CT is optional)
if exist(['./pCTMexFiles/' archdir],'dir')
    addpath(['./pCTMexFiles/' archdir]);
end
clear archdir


if ispc
    [user, sys]=memory;
    
    if sys.PhysicalMemory.Total<9000000000 % 8Gb
        warning('Your Computer has 8Gb or less of RAM memory. Using image sizes of higher than 512^3 is not recommended (most likely not possible)')
    end
    
    if sys.PhysicalMemory.Total<2500000000 % 2Gb
        warning('Your Computer has 2Gb or less of RAM memory. Using image sizes of higher than 256^3 is not recommended (most likely not possible)')
    end
else
    warning('TIGRE needs a big amount of memory, be careful when running big images.')
end

clear all;
