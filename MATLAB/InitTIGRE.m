function InitTIGRE(varargin)
%INITTIGRE Initialize TIGRE and select its native compute backend.
%   InitTIGRE() initializes the legacy CUDA backend.
%   InitTIGRE('backend',BACKEND) selects "cuda", "cpu", or "metal" for
%   subsequent TIGRE calls. The selection can be changed by calling
%   InitTIGRE again and survives CLEAR because it is stored in root
%   application data.

%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License: Open Source under BSD.
%--------------------------------------------------------------------------

parser = inputParser;
parser.FunctionName = 'InitTIGRE';
addParameter(parser, 'backend', 'cuda', @isSupportedBackend);
parse(parser, varargin{:});
backend = lower(char(string(parser.Results.backend)));

matlabRoot = fileparts(mfilename('fullpath'));
addToolboxPaths(matlabRoot);

backendFolder = getBackendFolder(matlabRoot, backend);
if ~isfolder(backendFolder)
    error('TIGRE:InitTIGRE:BackendNotBuilt', ...
        ['The TIGRE %s backend is not built for %s.\n' ...
         'Build it with Compile(''backend'',''%s'') and run InitTIGRE again.'], ...
        upper(backend), computer('arch'), backend);
end

addpath(backendFolder, '-begin');
if strcmp(backend, 'cuda')
    addpath(fullfile(matlabRoot, 'pCTMexFiles', legacyPlatformFolder()), '-begin');
end
rehash;
entryPoints = getBackendEntryPoints(backend);
for index = 1:numel(entryPoints)
    if exist(entryPoints{index}, 'file') ~= 3
        error('TIGRE:InitTIGRE:BackendNotBuilt', ...
            ['The TIGRE %s backend binary (%s.%s) is unavailable.\n' ...
             'Build it with Compile(''backend'',''%s'') and run InitTIGRE again.'], ...
            upper(backend), entryPoints{index}, mexext, backend);
    end
end

state = struct('backend', backend, ...
               'root', matlabRoot, ...
               'binaryFolder', backendFolder);
setappdata(groot, 'TIGREBackendState', state);
end

function tf = isSupportedBackend(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    any(strcmpi(char(string(value)), {'cuda', 'cpu', 'metal'}));
end

function addToolboxPaths(matlabRoot)
folders = {
    'Algorithms'
    'Utilities'
    fullfile('Utilities', 'Quality_measures')
    fullfile('Utilities', 'IO', 'VarianCBCT')
    fullfile('Utilities', 'IO', 'Nikon')
    fullfile('Utilities', 'IO', 'Dicom')
    fullfile('Utilities', 'IO', 'Bruker')
    fullfile('Utilities', 'IO', 'YXLON')
    fullfile('Utilities', 'IO', 'Diondo')
    fullfile('Utilities', 'IO', 'DXChange')
    fullfile('Utilities', 'GPU')
    fullfile('Utilities', 'DE')
    'Demos'
    'Colormaps'
    fullfile('Third_party_tools', 'arrow3d')
    fullfile('Third_party_tools', 'sec2hours')
    fullfile('Third_party_tools', 'readMHD')
    fullfile('Third_party_tools', 'Inpaint_nans')
    };

for index = 1:numel(folders)
    addpath(fullfile(matlabRoot, folders{index}));
end
addpath(genpath(fullfile(matlabRoot, 'Test_data')));
end

function folder = getBackendFolder(matlabRoot, backend)
switch backend
    case 'cuda'
        folder = fullfile(matlabRoot, 'Mex_files', legacyPlatformFolder());
    case 'cpu'
        folder = fullfile(matlabRoot, 'Mex_files', 'cpu', computer('arch'));
    case 'metal'
        if ~ismac || ~strcmp(computer('arch'), 'maca64')
            error('TIGRE:InitTIGRE:UnsupportedPlatform', ...
                'The Metal backend requires Apple Silicon MATLAB on macOS.');
        end
        folder = fullfile(matlabRoot, 'Mex_files', 'metal', 'maca64');
end
end

function names = getBackendEntryPoints(backend)
switch backend
    case 'cuda'
        names = {'Ax_mex', 'Atb_mex', 'minTV', 'AwminTV', ...
            'tvDenoise', 'minPICCS', 'AddNoise', ...
            'getGpuCount_mex', 'getGpuName_mex'};
    case 'cpu'
        names = {'tigre_cpu_mex'};
    case 'metal'
        names = {'tigre_metal_mex'};
end
end

function platform = legacyPlatformFolder()
if ispc
    if contains(computer('arch'), '64')
        platform = 'win64';
    else
        platform = 'win32';
    end
elseif ismac
    if contains(computer('arch'), '64')
        platform = 'mac64';
    else
        platform = 'mac32';
    end
else
    if contains(computer('arch'), '64')
        platform = 'linux64';
    else
        platform = 'linux32';
    end
end
end
