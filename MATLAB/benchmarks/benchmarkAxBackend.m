function result = benchmarkAxBackend(backend, varargin)
%BENCHMARKAXBACKEND Benchmark warmed Ax calls on one native backend.
%   RESULT = BENCHMARKAXBACKEND(BACKEND) benchmarks the Metal or CUDA
%   implementation using deterministic cone, parallel, helical, DBT, and
%   small-overhead cases.
%
%   RESULT = BENCHMARKAXBACKEND(..., 'OutputFile', FILE) saves RESULT to a
%   version 7.3 MAT file. Name-value options are:
%     Repetitions        Timed samples per projector and case (default 5).
%     Warmups            Untimed calls before each measurement (default 2).
%     OutputFile         MAT-file destination (default does not save).
%     RequestId          Identifier stored in the result metadata.
%     StoreOutputs       Retain the final output from every case (default true).
%     CompareCudaFixtures
%                       Record errors against existing fixtures (default true).
%     RunCgls            Time plot-free d13 CGLS (default false).
%     CglsIterations     Requested CGLS iteration count (default 20).

backend = validatestring(backend, {'metal', 'cuda'}, ...
    mfilename, 'backend', 1);

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'Repetitions', 5, @isPositiveInteger);
addParameter(parser, 'Warmups', 2, @isPositiveInteger);
addParameter(parser, 'OutputFile', "", @isTextScalar);
addParameter(parser, 'RequestId', "", @isTextScalar);
addParameter(parser, 'StoreOutputs', true, @isLogicalScalar);
addParameter(parser, 'CompareCudaFixtures', true, @isLogicalScalar);
addParameter(parser, 'RunCgls', false, @isLogicalScalar);
addParameter(parser, 'CglsIterations', 20, @isPositiveInteger);
parse(parser, varargin{:});
options = parser.Results;
if options.Repetitions < 5
    error('TIGRE:BenchmarkAx:TooFewRepetitions', ...
        'Repetitions must be at least five.');
end
if options.Warmups < 2
    error('TIGRE:BenchmarkAx:TooFewWarmups', ...
        'Warmups must be at least two.');
end
if ~options.StoreOutputs && ...
        (options.CompareCudaFixtures || options.RunCgls)
    error('TIGRE:BenchmarkAx:OutputsRequired', ...
        ['StoreOutputs must be true when comparing CUDA fixtures or ' ...
         'running CGLS.']);
end

matlabRoot = fileparts(fileparts(mfilename('fullpath')));
repositoryRoot = fileparts(matlabRoot);
addpath(matlabRoot);
InitTIGRE('backend', backend);

requestId = string(options.RequestId);
if strlength(requestId) == 0
    timestamp = datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyyMMdd''T''HHmmss''Z''');
    requestId = string(timestamp) + "-" + shortRevision(repositoryRoot);
end

previousFigureVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
figureCleanup = onCleanup( ...
    @() set(groot, 'DefaultFigureVisible', previousFigureVisibility));

result = struct;
result.schemaVersion = 1;
result.requestId = char(requestId);
result.backend = backend;
result.startedAtUtc = char(datetime('now', 'TimeZone', 'UTC'));
result.environment = environmentMetadata(repositoryRoot, backend);
result.configuration = struct( ...
    'repetitions', options.Repetitions, ...
    'warmups', options.Warmups, ...
    'storeOutputs', options.StoreOutputs, ...
    'compareCudaFixtures', options.CompareCudaFixtures, ...
    'runCgls', options.RunCgls, ...
    'cglsIterations', options.CglsIterations);
result.sourceFiles = sourceFileMetadata(repositoryRoot);

cases = benchmarkCases();
projectors = {'Siddon', 'interpolated'};
caseResults = repmat(emptyCaseResult(), size(cases));
fprintf('Ax benchmark %s on %s (%s)\n', requestId, upper(backend), ...
    result.environment.deviceNames{1});
fprintf('Warmups: %d, timed repetitions: %d\n', ...
    options.Warmups, options.Repetitions);

for caseIndex = 1:numel(cases)
    caseData = cases(caseIndex);
    caseResult = emptyCaseResult();
    caseResult.name = caseData.name;
    caseResult.imageSize = size3(caseData.image);
    caseResult.detectorSize = double(caseData.geometry.nDetector(:).');
    caseResult.angleSize = size(caseData.angles);
    caseResult.mode = char(string(caseData.geometry.mode));
    caseResult.projectors = repmat(emptyProjectorResult(), ...
        1, numel(projectors));

    for projectorIndex = 1:numel(projectors)
        projector = projectors{projectorIndex};
        output = [];
        for warmupIndex = 1:options.Warmups
            output = Ax( ...
                caseData.image, caseData.geometry, caseData.angles, ...
                projector);
            validateOutput(output, caseData);
        end

        times = zeros(1, options.Repetitions);
        for repetition = 1:options.Repetitions
            synchronizeBackend(backend);
            startTime = tic;
            output = Ax( ...
                caseData.image, caseData.geometry, caseData.angles, ...
                projector);
            synchronizeBackend(backend);
            times(repetition) = toc(startTime);
            validateOutput(output, caseData);
        end

        projectorResult = emptyProjectorResult();
        projectorResult.name = projector;
        projectorResult.timesSeconds = times;
        projectorResult.medianSeconds = median(times);
        projectorResult.minimumSeconds = min(times);
        projectorResult.maximumSeconds = max(times);
        projectorResult.medianAbsoluteDeviationSeconds = ...
            median(abs(times - projectorResult.medianSeconds));
        projectorResult.outputSize = size3(output);
        projectorResult.outputClass = class(output);
        projectorResult.allFinite = all(isfinite(output), 'all');
        projectorResult.outputSummary = arraySummary(output);
        if options.StoreOutputs
            projectorResult.output = output;
        end
        caseResult.projectors(projectorIndex) = projectorResult;

        fprintf('  %-22s %-12s median %.6f s', ...
            caseData.name, projector, projectorResult.medianSeconds);
        fprintf(' [%.6f, %.6f], MAD %.6f\n', ...
            projectorResult.minimumSeconds, ...
            projectorResult.maximumSeconds, ...
            projectorResult.medianAbsoluteDeviationSeconds);
    end
    caseResults(caseIndex) = caseResult;
end
result.cases = caseResults;

if options.CompareCudaFixtures
    result.cudaReferenceMetrics = compareCudaFixtures( ...
        cases, caseResults, matlabRoot);
else
    result.cudaReferenceMetrics = struct([]);
end

if options.RunCgls
    d13Index = find(strcmp({cases.name}, 'd13_helical'), 1);
    siddonIndex = find(strcmp( ...
        {caseResults(d13Index).projectors.name}, 'Siddon'), 1);
    projections = caseResults(d13Index).projectors(siddonIndex).output;
    cglsHandle = @() runCgls( ...
        projections, cases(d13Index).geometry, cases(d13Index).angles, ...
        options.CglsIterations);
    synchronizeBackend(backend);
    startTime = tic;
    [reconstruction, residual] = cglsHandle();
    synchronizeBackend(backend);
    result.cgls = struct( ...
        'requestedIterations', options.CglsIterations, ...
        'elapsedSeconds', toc(startTime), ...
        'residualLength', numel(residual), ...
        'reconstructionSize', size3(reconstruction), ...
        'reconstructionClass', class(reconstruction), ...
        'allFinite', all(isfinite(reconstruction), 'all') && ...
            all(isfinite(residual), 'all'), ...
        'residual', residual);
    if options.StoreOutputs
        result.cgls.reconstruction = reconstruction;
    end
    fprintf('  d13 CGLS(%d): %.6f s, residual length %d\n', ...
        options.CglsIterations, result.cgls.elapsedSeconds, ...
        result.cgls.residualLength);
else
    result.cgls = struct([]);
end

result.completedAtUtc = char(datetime('now', 'TimeZone', 'UTC'));
validateBenchmarkResult(result);

outputFile = string(options.OutputFile);
if strlength(outputFile) > 0
    outputFile = absolutePath(outputFile, repositoryRoot);
    outputFolder = fileparts(outputFile);
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    save(outputFile, 'result', '-v7.3');
    fprintf('Saved benchmark result: %s\n', outputFile);
end

clear figureCleanup
end

function cases = benchmarkCases()
cases = [ ...
    d13Case(), ...
    circularConeCase(), ...
    parallelCase(), ...
    staticDetectorCase(), ...
    smallOverheadCase()];
end

function value = d13Case()
geometry = struct;
geometry.DSD = 1536;
geometry.DSO = 1000;
geometry.nDetector = [512; 512];
geometry.dDetector = [0.8; 0.8];
geometry.sDetector = geometry.nDetector .* geometry.dDetector;
geometry.nVoxel = [128; 128; 512];
geometry.sVoxel = [256; 256; 1024];
geometry.dVoxel = geometry.sVoxel ./ geometry.nVoxel;
geometry.offDetector = [0; 0];
geometry.COR = 0;
geometry.mode = 'cone';
geometry.accuracy = 0.5;
angles = linspace(0, 2 * pi, 100);
angles = [angles angles angles];
geometry.offOrigin = zeros(3, numel(angles));
geometry.offOrigin(3, :) = linspace(-384, 384, numel(angles));
value = makeCase( ...
    'd13_helical', single(headPhantom(geometry.nVoxel)), ...
    geometry, angles);
end

function value = circularConeCase()
geometry = benchmarkGeometry([96; 96; 64], [160; 120], 'cone');
angles = linspace(0, 2 * pi, 91);
angles(end) = [];
geometry.offOrigin = [ ...
    1.2 * sin(angles); ...
    0.8 * cos(angles); ...
    linspace(-3, 3, numel(angles))];
geometry.offDetector = [1.5 * sin(angles); 0.9 * cos(angles)];
geometry.COR = 0.3 * sin(angles);
geometry.rotDetector = [ ...
    0.03 * sin(angles); ...
    0.02 * cos(angles); ...
    0.025 * sin(2 * angles)];
value = makeCase( ...
    'circular_cone', benchmarkPhantom(geometry.nVoxel), geometry, angles);
end

function value = parallelCase()
geometry = benchmarkGeometry([96; 96; 64], [160; 96], 'parallel');
angles = linspace(0, pi, 91);
angles(end) = [];
geometry.offOrigin = [ ...
    0.8 * sin(angles); ...
    0.5 * cos(angles); ...
    linspace(-2, 2, numel(angles))];
geometry.offDetector = [0.7 * sin(angles); 0.4 * cos(angles)];
geometry.COR = 0.2 * sin(angles);
value = makeCase( ...
    'parallel', benchmarkPhantom(geometry.nVoxel), geometry, angles);
end

function value = staticDetectorCase()
geometry = benchmarkGeometry([64; 128; 48], [160; 96], 'cone');
geometry.DSO = 620;
geometry.DSD = 660;
geometry.sVoxel = [32; 64; 24];
geometry.dVoxel = geometry.sVoxel ./ geometry.nVoxel;
geometry.dDetector = [0.4; 0.4];
geometry.sDetector = geometry.nDetector .* geometry.dDetector;
geometry.offOrigin = [12; 0; 12];
geometry.offDetector = [0; geometry.sDetector(2) / 2];
angles = linspace(-pi / 12, pi / 12, 31);
geometry = staticDetectorGeo(geometry, angles);
value = makeCase( ...
    'static_detector_dbt', benchmarkPhantom(geometry.nVoxel), ...
    geometry, angles);
end

function value = smallOverheadCase()
geometry = benchmarkGeometry([8; 8; 8], [8; 8], 'cone');
angles = [0, pi / 3, 2 * pi / 3];
value = makeCase( ...
    'small_overhead', benchmarkPhantom(geometry.nVoxel), ...
    geometry, angles);
end

function geometry = benchmarkGeometry(nVoxel, nDetector, mode)
geometry = defaultGeometry( ...
    'nVoxel', nVoxel, 'nDetector', nDetector, 'mode', mode);
geometry.nVoxel = double(nVoxel);
geometry.nDetector = double(nDetector);
geometry.sVoxel = double(nVoxel);
geometry.dVoxel = geometry.sVoxel ./ geometry.nVoxel;
geometry.dDetector = [1; 1];
geometry.sDetector = geometry.nDetector .* geometry.dDetector;
geometry.DSO = 400;
geometry.DSD = 600;
geometry.accuracy = 0.5;
end

function image = benchmarkPhantom(dimensions)
[x, y, z] = ndgrid( ...
    linspace(-1, 1, dimensions(1)), ...
    linspace(-1, 1, dimensions(2)), ...
    linspace(-1, 1, dimensions(3)));
image = single( ...
    0.8 * (x.^2 + 1.2 * y.^2 + 0.8 * z.^2 < 0.58) + ...
    0.2 * ((x + 0.35).^2 + (y - 0.2).^2 + z.^2 < 0.08));
end

function value = makeCase(name, image, geometry, angles)
value = struct( ...
    'name', name, ...
    'image', image, ...
    'geometry', geometry, ...
    'angles', angles);
end

function metrics = compareCudaFixtures(cases, caseResults, matlabRoot)
metrics = struct( ...
    'caseName', {}, ...
    'projector', {}, ...
    'relativeL2Error', {}, ...
    'maximumAbsoluteError', {}, ...
    'rmse', {}, ...
    'referenceFinite', {}, ...
    'actualFinite', {});

d13FixtureFile = fullfile( ...
    matlabRoot, 'tests', 'fixtures', 'd13_cuda_reference.mat');
if isfile(d13FixtureFile)
    loaded = load(d13FixtureFile, 'trial');
    d13CaseIndex = find(strcmp({cases.name}, 'd13_helical'), 1);
    siddonIndex = find(strcmp( ...
        {caseResults(d13CaseIndex).projectors.name}, 'Siddon'), 1);
    actual = caseResults(d13CaseIndex).projectors(siddonIndex).output;
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'd13_helical', 'Siddon', loaded.trial.cuda.projections, actual);
    clear loaded actual
end

smallFixtureFile = fullfile( ...
    matlabRoot, 'tests', 'fixtures', 'cuda_reference.mat');
if isfile(smallFixtureFile)
    loaded = load(smallFixtureFile, 'fixture');
    InitTIGRE('backend', char(tigreBackend()));
    actual = Ax( ...
        loaded.fixture.input.coneImage, ...
        loaded.fixture.input.coneGeometry, ...
        loaded.fixture.input.coneAngles, 'Siddon');
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'cuda_fixture_cone', 'Siddon', ...
        loaded.fixture.cuda.coneSiddon, actual);
    actual = Ax( ...
        loaded.fixture.input.coneImage, ...
        loaded.fixture.input.coneGeometry, ...
        loaded.fixture.input.coneAngles, 'interpolated');
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'cuda_fixture_cone', 'interpolated', ...
        loaded.fixture.cuda.coneInterpolated, actual);
    actual = Ax( ...
        loaded.fixture.input.parallelImage, ...
        loaded.fixture.input.parallelGeometry, ...
        loaded.fixture.input.parallelAngles, 'Siddon');
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'cuda_fixture_parallel', 'Siddon', ...
        loaded.fixture.cuda.parallelSiddon, actual);
    actual = Ax( ...
        loaded.fixture.input.parallelImage, ...
        loaded.fixture.input.parallelGeometry, ...
        loaded.fixture.input.parallelAngles, 'interpolated');
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'cuda_fixture_parallel', 'interpolated', ...
        loaded.fixture.cuda.parallelInterpolated, actual);
    clear loaded actual
end

geometryFixtureFile = fullfile( ...
    matlabRoot, 'tests', 'fixtures', 'cuda_geometry_reference.mat');
if isfile(geometryFixtureFile)
    loaded = load(geometryFixtureFile, 'fixture');
    caseData = loaded.fixture.cases.staticDetectorDbt;
    actual = Ax( ...
        caseData.image, caseData.geometry, caseData.angles, 'Siddon');
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'cuda_fixture_static_detector_dbt', 'Siddon', ...
        caseData.cuda.siddon, actual);
    actual = Ax( ...
        caseData.image, caseData.geometry, caseData.angles, 'interpolated');
    metrics(end + 1) = errorMetrics( ... %#ok<AGROW>
        'cuda_fixture_static_detector_dbt', 'interpolated', ...
        caseData.cuda.interpolated, actual);
end

for metricIndex = 1:numel(metrics)
    fprintf('  CUDA error %-33s %-12s relL2 %.3e, max %.3e, RMSE %.3e\n', ...
        metrics(metricIndex).caseName, metrics(metricIndex).projector, ...
        metrics(metricIndex).relativeL2Error, ...
        metrics(metricIndex).maximumAbsoluteError, ...
        metrics(metricIndex).rmse);
end
end

function value = errorMetrics(caseName, projector, reference, actual)
if ~isequal(size(reference), size(actual))
    error('TIGRE:BenchmarkAx:ReferenceSizeMismatch', ...
        'CUDA reference size differs for case %s, projector %s.', ...
        caseName, projector);
end

chunkSize = 4 * 1024 * 1024;
sumSquaredError = 0;
sumSquaredReference = 0;
maximumAbsoluteError = 0;
referenceFinite = true;
actualFinite = true;
elementCount = numel(reference);
for first = 1:chunkSize:elementCount
    last = min(elementCount, first + chunkSize - 1);
    referenceChunk = double(reference(first:last));
    actualChunk = double(actual(first:last));
    difference = actualChunk - referenceChunk;
    sumSquaredError = sumSquaredError + sum(difference .* difference);
    sumSquaredReference = sumSquaredReference + ...
        sum(referenceChunk .* referenceChunk);
    maximumAbsoluteError = max( ...
        maximumAbsoluteError, max(abs(difference), [], 'all'));
    referenceFinite = referenceFinite && all(isfinite(referenceChunk));
    actualFinite = actualFinite && all(isfinite(actualChunk));
end

value = struct( ...
    'caseName', caseName, ...
    'projector', projector, ...
    'relativeL2Error', sqrt(sumSquaredError) / ...
        max(sqrt(sumSquaredReference), eps), ...
    'maximumAbsoluteError', maximumAbsoluteError, ...
    'rmse', sqrt(sumSquaredError / max(elementCount, 1)), ...
    'referenceFinite', referenceFinite, ...
    'actualFinite', actualFinite);
end

function metadata = environmentMetadata(repositoryRoot, backend)
metadata = struct( ...
    'matlabVersion', version, ...
    'computer', computer, ...
    'operatingSystem', systemOutput('sw_vers'), ...
    'deviceNames', {getGpuNames()}, ...
    'sourceVersion', sourceVersion(repositoryRoot), ...
    'sourceStatus', sourceStatus(repositoryRoot), ...
    'backendBinary', which(backendBinaryName(backend)));
end

function files = sourceFileMetadata(repositoryRoot)
relativePaths = { ...
    fullfile('MATLAB', 'benchmarks', 'benchmarkAxBackend.m'), ...
    fullfile('MATLAB', 'Utilities', 'Ax.m'), ...
    fullfile('MATLAB', 'Utilities', 'checkGeo.m'), ...
    fullfile('MATLAB', 'Algorithms', 'CGLS.m'), ...
    fullfile('Common', 'Geometry', 'tigre_geometry.hpp'), ...
    fullfile('Common', 'Metal', 'tigre_metal_mex.mm'), ...
    fullfile('Common', 'Metal', 'tigre_kernels.metal')};
files = repmat(struct('path', '', 'sha256', ''), size(relativePaths));
for fileIndex = 1:numel(relativePaths)
    files(fileIndex).path = relativePaths{fileIndex};
    files(fileIndex).sha256 = fileSha256( ...
        fullfile(repositoryRoot, relativePaths{fileIndex}));
end
end

function value = fileSha256(filename)
if ispc
    escaped = strrep(filename, '''', '''''');
    command = sprintf([ ...
        'powershell -NoProfile -Command ' ...
        '"(Get-FileHash -Algorithm SHA256 -LiteralPath ''%s'').Hash"'], ...
        escaped);
else
    command = sprintf('shasum -a 256 "%s"', ...
        strrep(filename, '"', '\"'));
end
[status, output] = system(command);
if status ~= 0
    value = 'unavailable';
    return
end
tokens = regexp(strtrim(output), '^[0-9A-Fa-f]+', 'match', 'once');
if isempty(tokens)
    value = 'unavailable';
else
    value = lower(tokens);
end
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

function value = shortRevision(repositoryRoot)
revision = sourceVersion(repositoryRoot);
value = string(revision(1:min(7, strlength(revision))));
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

function value = systemOutput(command)
[status, output] = system(command);
if status == 0
    value = strtrim(output);
else
    value = 'unavailable';
end
end

function name = backendBinaryName(backend)
if strcmp(backend, 'cuda')
    name = 'Ax_mex';
else
    name = 'tigre_metal_mex';
end
end

function synchronizeBackend(backend)
if strcmp(backend, 'cuda')
    wait(gpuDevice);
end
end

function [reconstruction, residual] = runCgls( ...
        projections, geometry, angles, iterations)
[reconstruction, residual] = CGLS( ...
    projections, geometry, angles, iterations, 'Verbose', false);
end

function validateOutput(output, caseData)
expectedSize = [ ...
    caseData.geometry.nDetector(2), ...
    caseData.geometry.nDetector(1), ...
    size(caseData.angles, 2)];
if ~isa(output, 'single') || ~isequal(size3(output), expectedSize) || ...
        ~all(isfinite(output), 'all')
    error('TIGRE:BenchmarkAx:InvalidOutput', ...
        'Ax returned an invalid output for case %s.', caseData.name);
end
end

function validateBenchmarkResult(result)
assert(result.schemaVersion == 1);
assert(~isempty(result.requestId));
assert(any(strcmp(result.backend, {'metal', 'cuda'})));
assert(numel(result.cases) == 5);
for caseIndex = 1:numel(result.cases)
    assert(numel(result.cases(caseIndex).projectors) == 2);
    for projectorIndex = 1:2
        projector = result.cases(caseIndex).projectors(projectorIndex);
        assert(numel(projector.timesSeconds) >= 5);
        assert(all(isfinite(projector.timesSeconds)));
        assert(all(projector.timesSeconds > 0));
        assert(projector.allFinite);
    end
end
if ~isempty(result.cgls)
    assert(result.cgls.allFinite);
    assert(result.cgls.residualLength >= ...
        result.cgls.requestedIterations);
end
end

function value = arraySummary(array)
values = double(array(:));
value = struct( ...
    'minimum', min(values), ...
    'maximum', max(values), ...
    'mean', mean(values), ...
    'l2Norm', norm(values));
end

function value = emptyCaseResult()
value = struct( ...
    'name', '', ...
    'imageSize', [], ...
    'detectorSize', [], ...
    'angleSize', [], ...
    'mode', '', ...
    'projectors', struct([]));
end

function value = emptyProjectorResult()
value = struct( ...
    'name', '', ...
    'timesSeconds', [], ...
    'medianSeconds', [], ...
    'minimumSeconds', [], ...
    'maximumSeconds', [], ...
    'medianAbsoluteDeviationSeconds', [], ...
    'outputSize', [], ...
    'outputClass', '', ...
    'allFinite', false, ...
    'outputSummary', struct, ...
    'output', []);
end

function value = size3(array)
value = size(array);
value(end + 1:3) = 1;
end

function value = absolutePath(filename, baseFolder)
filename = char(filename);
if ispc
    isAbsolute = ~isempty(regexp(filename, ...
        '^[A-Za-z]:[\\/]|^\\\\', 'once'));
else
    isAbsolute = startsWith(filename, filesep);
end
if isAbsolute
    value = filename;
else
    value = fullfile(baseFolder, filename);
end
end

function value = isPositiveInteger(input)
value = isnumeric(input) && isscalar(input) && isfinite(input) && ...
    input > 0 && input == fix(input);
end

function value = isTextScalar(input)
value = (ischar(input) && isrow(input)) || ...
    (isstring(input) && isscalar(input));
end

function value = isLogicalScalar(input)
value = islogical(input) && isscalar(input);
end
