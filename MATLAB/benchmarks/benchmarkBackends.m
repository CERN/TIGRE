function results = benchmarkBackends()
%BENCHMARKBACKENDS Benchmark warmed CPU and Metal MVP native operations.
%   RESULTS = benchmarkBackends() returns median-style TIMEIT measurements.
%   This benchmark is informational and deliberately has no speed gate.

matlabRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(matlabRoot);

geometry = defaultGeometry( ...
    'nVoxel', [24; 24; 16], 'nDetector', [32; 20], 'mode', 'cone');
geometry.sVoxel = [24; 24; 16];
geometry.dVoxel = geometry.sVoxel ./ geometry.nVoxel;
geometry.dDetector = [1; 1];
geometry.sDetector = geometry.nDetector .* geometry.dDetector;
geometry.DSO = 60;
geometry.DSD = 90;
geometry.accuracy = 0.75;
angles = linspace(0, 2 * pi, 17);
angles(end) = [];
[x, y, z] = ndgrid( ...
    linspace(-1, 1, geometry.nVoxel(1)), ...
    linspace(-1, 1, geometry.nVoxel(2)), ...
    linspace(-1, 1, geometry.nVoxel(3)));
image = single((x.^2 + y.^2 + z.^2) < 0.6);

backends = ["cpu"; "metal"];
axSeconds = zeros(size(backends));
atbSeconds = zeros(size(backends));
tvSeconds = zeros(size(backends));
for index = 1:numel(backends)
    InitTIGRE('backend', backends(index));
    gpuIds = GpuIds();
    projections = Ax(image, geometry, angles, 'interpolated');
    Atb(projections, geometry, angles, 'matched');
    minimizeTV(image, 0.05, 2, 'gpuids', gpuIds);

    axSeconds(index) = timeit( ...
        @() Ax(image, geometry, angles, 'interpolated'));
    atbSeconds(index) = timeit( ...
        @() Atb(projections, geometry, angles, 'matched'));
    tvSeconds(index) = timeit( ...
        @() minimizeTV(image, 0.05, 2, 'gpuids', gpuIds));
end

results = table(backends, axSeconds, atbSeconds, tvSeconds, ...
    'VariableNames', {'Backend', 'AxSeconds', 'AtbSeconds', 'TvSeconds'});
disp(results);
end
