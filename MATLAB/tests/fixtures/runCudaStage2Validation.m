function [fixtureFile, logFile] = runCudaStage2Validation()
%RUNCUDASTAGE2VALIDATION Build CUDA, create fixtures, and run Stage 2 demos.

fixtureFolder = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(fileparts(fixtureFolder));
repositoryRoot = fileparts(matlabRoot);
fixtureFile = fullfile(fixtureFolder, 'cuda_geometry_reference.mat');
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
logFile = fullfile( ...
    fixtureFolder, ['cuda_stage2_validation_', timestamp, '.log']);

previousVisibility = get(groot, 'DefaultFigureVisible');
previousFigures = findall(groot, 'Type', 'figure');
previousPath = path;
previousFolder = pwd;
previousRng = rng;
hadBackend = isappdata(groot, 'TIGREBackendState');
previousBackend = getappdata(groot, 'TIGREBackendState');
cleanup = onCleanup(@() finishValidation(previousVisibility, previousFigures, ...
    previousPath, previousFolder, previousRng, hadBackend, previousBackend));
set(groot, 'DefaultFigureVisible', 'off');
diary(logFile);

try
    addpath(matlabRoot);
    addpath(fixtureFolder);
    addpath(fullfile(matlabRoot, 'tests', 'backend'));

    fprintf('TIGRE Stage 2 CUDA validation\n');
    fprintf('Started: %s\n', char(datetime('now')));
    fprintf('Repository: %s\n', repositoryRoot);
    fprintf('Source version: %s\n', sourceVersion(repositoryRoot));
    fprintf('Source status:\n%s\n', sourceStatus(repositoryRoot));
    fprintf('MATLAB: %s\n', version);
    fprintf('Computer: %s\n', computer);

    clear Ax_mex Atb_mex minTV AwminTV tvDenoise minPICCS AddNoise
    clear getGpuCount_mex getGpuName_mex
    Compile('backend', 'cuda');
    InitTIGRE('backend', 'cuda');
    fprintf('CUDA devices:\n');
    disp(getGpuNames());

    generateCudaGeometryReferenceFixtures(fixtureFile);
    loaded = load(fixtureFile, 'fixture');
    validateCudaGeometryReferenceFixture(loaded.fixture);
    fprintf('Fixture validation passed: %s\n', fixtureFile);

    demos = { ...
        'd13_HelicalGeometry.m', ...
        'd14_Offsets.m', ...
        'd17_DetectorRotation.m', ...
        'd18_ArbitraryAxisOfRotation.m', ...
        'd19_DBTexample.m', ...
        'd22_ListGPUs.m'};
    for demoIndex = 1:numel(demos)
        runDemo(demos{demoIndex});
    end

    fprintf('Completed: %s\n', char(datetime('now')));
    fprintf('Stage 2 CUDA validation passed.\n');
catch exception
    fprintf(2, '%s\n', getReport( ...
        exception, 'extended', 'hyperlinks', 'off'));
    rethrow(exception);
end

clear cleanup
end

function runDemo(demoName)
fprintf('Running %s\n', demoName);
startTime = tic;
backendtest.runIsolated('demo', 'cuda', demoName);
fprintf('Passed %s in %.3f seconds\n', demoName, toc(startTime));
end

function value = sourceVersion(repositoryRoot)
[status, output] = system(sprintf( ...
    'git -C "%s" rev-parse HEAD', repositoryRoot));
if status == 0
    value = strtrim(output);
else
    value = 'unavailable';
end
end

function value = sourceStatus(repositoryRoot)
[status, output] = system(sprintf( ...
    'git -C "%s" status --short', repositoryRoot));
if status ~= 0
    value = 'unavailable';
elseif isempty(strtrim(output))
    value = 'clean';
else
    value = strtrim(output);
end
end

function finishValidation(previousVisibility, previousFigures, ...
        previousPath, previousFolder, previousRng, hadBackend, previousBackend)
delete(setdiff(findall(groot, 'Type', 'figure'), previousFigures));
set(groot, 'DefaultFigureVisible', previousVisibility);
rng(previousRng);
cd(previousFolder);
path(previousPath);
if hadBackend
    setappdata(groot, 'TIGREBackendState', previousBackend);
elseif isappdata(groot, 'TIGREBackendState')
    rmappdata(groot, 'TIGREBackendState');
end
diary off
end
