function [backend, reason]=select_backend(requested)
% SELECT_BACKEND  Choose the GPU backend that TIGRE is compiled for.
%
%   [backend, reason]=select_backend(requested) returns the backend ('cuda',
%   'hip', 'sycl' or 'metal') and why it was chosen. requested is the name
%   given to Compile -backend, or empty. If it is empty, the TIGRE_BACKEND
%   environment variable is used. If that is not set either, the first
%   backend whose compiler is installed is used, in the order above. This is
%   the same rule as the Python build (setup.py).
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
backends={'cuda','hip','sycl','metal'};
if nargin<1
    requested='';
end
requested=lower(strtrim(char(requested)));
reason='set by Compile -backend';
if isempty(requested)
    requested=lower(strtrim(getenv('TIGRE_BACKEND')));
    reason='set by TIGRE_BACKEND';
end
if ~isempty(requested)
    if ~any(strcmp(requested,backends))
        error('TIGRE:Compile:backend','Backend "%s" is not known. Use one of: %s.',requested,strjoin(backends,', '));
    end
    backend=requested;
    return
end
detectors={@detect_cuda,@detect_hip,@detect_sycl,@detect_metal};
for i=1:numel(backends)
    reason=detectors{i}();
    if ~isempty(reason)
        backend=backends{i};
        return
    end
end
error('TIGRE:Compile:noCompiler',['TIGRE could not find a GPU compiler. Install the CUDA toolkit (NVIDIA GPUs) ' ...
    'and make sure nvcc is on the PATH, or tell MATLAB where it is with setenv(''CUDA_PATH'',''your path''). ' ...
    'To choose a backend yourself, run Compile -backend cuda (one of: %s).'],strjoin(backends,', '));
end

function reason=detect_cuda()
reason='';
if ~isempty(getenv('CUDA_PATH'))
    reason='CUDA_PATH is set';
elseif ~isempty(getenv('CUDA_HOME'))
    reason='CUDA_HOME is set';
elseif on_path('nvcc')
    reason='nvcc is on the PATH';
elseif ispc && exist('C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/','dir')
    reason='CUDA is installed in C:/Program Files/NVIDIA GPU Computing Toolkit';
elseif ~ispc && exist('/usr/local/cuda','dir')
    reason='CUDA is installed in /usr/local/cuda';
end
end

function reason=detect_hip()
reason='';
if ~isempty(getenv('ROCM_PATH'))
    reason='ROCM_PATH is set';
elseif ~isempty(getenv('HIP_PATH'))
    reason='HIP_PATH is set';
elseif on_path('hipcc')
    reason='hipcc is on the PATH';
elseif ~ispc && exist('/opt/rocm/bin/hipcc','file')
    reason='ROCm is installed in /opt/rocm';
end
end

function reason=detect_sycl()
reason='';
if ispc
    compiler='icx';
else
    compiler='icpx';
end
if on_path(compiler)
    reason=[compiler,' (oneAPI DPC++) is on the PATH'];
end
end

function reason=detect_metal()
reason='';
if ismac && strcmp(computer('arch'),'maca64')
    reason='this is a Mac with Apple Silicon';
end
end

function found=on_path(program)
if ispc
    [status,~]=system(['where ',program]);
else
    [status,~]=system(['which ',program]);
end
found=(status==0);
end
