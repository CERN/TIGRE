function worker(requestFile)
%WORKER Entry point for isolated backend checks; never called in the desktop.

assert(batchStartupOptionUsed, 'TIGRE:Tests:BatchRequired', ...
    'This worker must run in a separate MATLAB -batch process.');
loaded = load(requestFile, 'request');
request = loaded.request;
matlabRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(matlabRoot);
InitTIGRE('backend', request.backend);
rng(42, 'twister');
set(groot, 'DefaultFigureVisible', 'off');
result = struct('passed', false, 'backend', request.backend);
try
    switch request.operation
        case 'clear'
            clearLocalWorkspace();
        case 'demo'
            result.arrays = runDemo(matlabRoot, request.detail);
        case 'recovery'
            checkRecovery(matlabRoot, request.folder, request.detail);
        case 'fail'
            evalin('base', 'clear');
            close('all', 'force');
            error('TIGRE:Tests:IntentionalChildFailure', 'Intentional isolation probe.');
        case 'noResult'
            return
    end
    assert(tigreBackend() == string(request.backend), ...
        'TIGRE:Tests:BackendChanged', 'The child changed backend.');
    result.passed = true;
    save(request.resultFile, 'result');
catch exception
    result.error = getReport(exception, 'extended', 'hyperlinks', 'off');
    save(request.resultFile, 'result');
    rethrow(exception);
end
end

function clearLocalWorkspace()
clear
% This intentionally clears only the disposable child process's base workspace.
evalin('base', 'clear');
close('all', 'force');
assert(tigreBackend() == "metal", ...
    'TIGRE:Tests:BackendLostAfterClear', 'Backend selection did not survive clear.');
end

function result = runDemo(matlabRoot, demoName)
demoFile = fullfile(matlabRoot, 'Demos', demoName);
switch demoName
    case 'd13_HelicalGeometry.m'
        projectionNames = {'data'};
        reconstructionNames = {'OSSARTimg', 'CGLSimg'};
    case 'd14_Offsets.m'
        projectionNames = {'projections', 'projections2', 'projections3'};
        reconstructionNames = {'res'};
    case 'd17_DetectorRotation.m'
        projectionNames = {'projections', 'projections2'};
        reconstructionNames = {'imgRotDet', 'imgnoRot'};
    case 'd18_ArbitraryAxisOfRotation.m'
        projectionNames = {'projections'};
        reconstructionNames = {'imgSIRT', 'imgCGLS'};
    case 'd19_DBTexample.m'
        projectionNames = {'projections'};
        reconstructionNames = {'recFDK', 'recSART'};
    case 'd22_ListGPUs.m'
        projectionNames = {};
        reconstructionNames = {};
    otherwise
        error('TIGRE:Tests:UnknownDemo', 'Unknown demo: %s', demoName);
end
% The disposable child's base workspace gives unchanged scripts their usual
% semantics. Explicit reads avoid collisions with functions such as HEAD.
evalin('base', sprintf('run(''%s'')', strrep(demoFile, '''', '''''')));
if strcmp(demoName, 'd22_ListGPUs.m')
    result = devicesDemo();
else
    projections = cellfun(@(name) evalin('base', name), ...
        projectionNames, 'UniformOutput', false);
    reconstructions = cellfun(@(name) evalin('base', name), ...
        reconstructionNames, 'UniformOutput', false);
    result = checkArrays(evalin('base', 'geo'), evalin('base', 'angles'), ...
        evalin('base', 'head'), projections, reconstructions);
end
end

function result = devicesDemo()
listGpuNames = evalin('base', 'listGpuNames');
gpuids = evalin('base', 'gpuids');
selectedGpuIds = evalin('base', 'selectedGpuIds');
assert(~isempty(listGpuNames) && all(cellfun(@ischar, listGpuNames)));
if tigreBackend() ~= "cuda"
    assert(isequal(gpuids.devices, int32(0)));
    assert(isequal(selectedGpuIds.devices, int32(0)));
end
result = struct('deviceCount', numel(listGpuNames));
end

function result = checkArrays(geometry, angles, phantom, projections, reconstructions)
imageSize = double(geometry.nVoxel(:).');
projectionSize = [double(geometry.nDetector([2, 1]).'), size(angles, 2)];
checkArray(phantom, imageSize);
for index = 1:numel(projections)
    checkArray(projections{index}, projectionSize);
end
for index = 1:numel(reconstructions)
    checkArray(reconstructions{index}, imageSize);
end
result = struct('imageSize', imageSize, 'projectionSize', projectionSize, ...
    'projectionCount', numel(projections), ...
    'reconstructionCount', numel(reconstructions));
end

function checkArray(value, expectedSize)
assert(isa(value, 'single') && isreal(value) && ~isempty(value), ...
    'TIGRE:Tests:InvalidDemoArray', 'Expected a nonempty real single array.');
assert(isequal([size(value, 1), size(value, 2), size(value, 3)], expectedSize), ...
    'TIGRE:Tests:InvalidDemoShape', 'Unexpected demo output dimensions.');
assert(all(isfinite(value), 'all'), ...
    'TIGRE:Tests:NonfiniteDemoArray', 'Demo output contains nonfinite values.');
end

function checkRecovery(matlabRoot, folder, failureMode)
originalBinary = which('tigre_metal_mex');
originalLibrary = fullfile(fileparts(originalBinary), 'tigre_kernels.metallib');
copyfile(originalBinary, folder);
libraryCopy = fullfile(folder, 'tigre_kernels.metallib');
if strcmp(failureMode, 'corrupt')
    % A MAT file is deliberately not a Metal shader library.
    copyfile(fullfile(folder, 'request.mat'), libraryCopy);
end
addpath(folder, '-begin');
rehash;
assert(strcmp(fileparts(which('tigre_metal_mex')), folder));
try
    getGpuNames();
    error('TIGRE:Tests:ExpectedLibraryFailure', 'Missing library error was not raised.');
catch exception
    assert(strcmp(exception.identifier, 'TIGRE:Metal:LibraryNotFound'), ...
        'Unexpected initialisation error: %s', exception.identifier);
end
copyfile(originalLibrary, libraryCopy);
names = getGpuNames();
assert(~isempty(names));
% No clear mex or backend reinitialisation between failure and recovery.
loaded = load(fullfile(matlabRoot, 'tests', 'fixtures', 'cuda_reference.mat'), 'fixture');
reference = loaded.fixture;
actual = Ax(reference.input.coneImage, reference.input.coneGeometry, ...
    reference.input.coneAngles, 'Siddon');
relativeError = norm(double(actual(:)) - double(reference.cuda.coneSiddon(:))) / ...
    norm(double(reference.cuda.coneSiddon(:)));
assert(relativeError <= 6e-3, 'TIGRE:Tests:RecoveryNumerics', ...
    'Recovered backend did not reproduce CUDA projections.');
end
