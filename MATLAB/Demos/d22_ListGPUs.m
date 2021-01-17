% List the names of installed GPUs
listGpuNames = getGpuNames();
disp('listGpuNames');
disp(listGpuNames);

% Choose one of them
%targetGpuName='';
targetGpuName = 'GeForce GTX 1060 6GB';
%targetGpuName = 'GeForce RTX 2080 Ti';
%targetGpuName = 'GeForce GTX 1070';

% You can get the list of GPU IDs
gpuids = getGpuIds(targetGpuName);
fprintf('Result of (%s)\n', targetGpuName);
fprintf('\t Number of devices: %d\n', length(gpuids));
disp(gpuids);
