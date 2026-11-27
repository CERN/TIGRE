function args = gpuarray_mex_flags()
% Extra mex arguments to build Ax and Atb with gpuArray support, or {} if
% the Parallel Computing Toolbox (its C API for gpuArrays, mxGPUArray) is not
% installed. Without it, TIGRE compiles as before and Ax/Atb run gpuArray
% inputs through a host copy, with a warning.
%
% The toolbox's header is in toolbox/parallel/gpu/extern/include (R2019b and
% later) or toolbox/distcomp/gpu/extern/include (older releases). The library
% is gpu.lib on Windows and libmwgpu on Linux, as mexcuda links them.
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
%--------------------------------------------------------------------------
args = {};
include_dir = '';
candidates = {fullfile(matlabroot,'toolbox','parallel','gpu','extern','include'), ...
              fullfile(matlabroot,'toolbox','distcomp','gpu','extern','include')};
for ii = 1:numel(candidates)
    if exist(fullfile(candidates{ii},'gpu','mxGPUArray.h'),'file')
        include_dir = candidates{ii};
        break
    end
end
if ispc
    lib_dir = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft');
    lib_name = 'gpu';
    lib_found = exist(fullfile(lib_dir,'gpu.lib'),'file') == 2;
else
    lib_dir = fullfile(matlabroot,'bin',computer('arch'));
    lib_name = 'mwgpu';
    lib_found = ~isempty(dir(fullfile(lib_dir,'libmwgpu.*')));
end
if isempty(include_dir) || ~lib_found
    disp('Parallel Computing Toolbox not found: compiling without gpuArray support in Ax and Atb.');
    return
end
disp('Parallel Computing Toolbox found: compiling Ax and Atb with gpuArray support.');
args = {'-DTIGRE_GPUARRAY', ['-I' include_dir], ['-L' lib_dir], ['-l' lib_name]};
end
