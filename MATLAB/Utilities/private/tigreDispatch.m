function varargout = tigreDispatch(operation, varargin)
%TIGREDISPATCH Route a native operation to the active TIGRE backend.

backend = char(tigreBackend());
switch lower(operation)
    case 'ax'
        validateDeviceSelection(backend, varargin{5});
        functionName = backendFunction(backend, 'Ax_mex');
        args = varargin(1:4);
    case 'atb'
        validateDeviceSelection(backend, varargin{5});
        functionName = backendFunction(backend, 'Atb_mex');
        args = varargin(1:4);
    case 'mintv'
        validateDeviceSelection(backend, varargin{4});
        functionName = backendFunction(backend, 'minTV');
        args = varargin(1:3);
    case 'awmintv'
        validateDeviceSelection(backend, varargin{5});
        functionName = backendFunction(backend, 'AwminTV');
        args = varargin(1:4);
    case 'tvdenoise'
        validateDeviceSelection(backend, varargin{4});
        functionName = backendFunction(backend, 'tvDenoise');
        args = varargin(1:3);
    case 'minpiccs'
        validateDeviceSelection(backend, varargin{6});
        functionName = backendFunction(backend, 'minPICCS');
        args = varargin(1:5);
    otherwise
        error('TIGRE:Backend:UnknownOperation', ...
            'Unknown TIGRE native operation "%s".', operation);
end

if exist(functionName, 'file') ~= 3
    error('TIGRE:Backend:BinaryUnavailable', ...
        ['The %s backend binary "%s" is not available. ' ...
         'Run Compile(''backend'',''%s'') and InitTIGRE(''backend'',''%s'').'], ...
        upper(backend), functionName, backend, backend);
end

if strcmp(backend, 'cuda')
    [varargout{1:nargout}] = feval(functionName, varargin{:});
else
    [varargout{1:nargout}] = feval(functionName, operation, args{:});
end
end

function functionName = backendFunction(backend, cudaFunction)
switch backend
    case 'cuda'
        functionName = cudaFunction;
    case 'cpu'
        functionName = 'tigre_cpu_mex';
    case 'metal'
        functionName = 'tigre_metal_mex';
    otherwise
        error('TIGRE:Backend:InvalidState', ...
            'Invalid active TIGRE backend "%s".', backend);
end
end

function validateDeviceSelection(backend, devices)
if strcmp(backend, 'cuda')
    return
end
if ~isa(devices, 'int32') || ~isequal(devices, int32(0))
    error('TIGRE:Backend:InvalidDevice', ...
        'The %s backend exposes one logical device with ID 0.', upper(backend));
end
end
