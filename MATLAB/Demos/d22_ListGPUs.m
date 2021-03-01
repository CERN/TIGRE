% List the names of installed GPUs
listGpuNames = getGpuNames();
disp('listGpuNames');
disp(listGpuNames);

% Choose one of them
%targetGpuName='';
%targetGpuName = 'GeForce GTX 1060 6GB';
targetGpuName = 'GeForce RTX 2080 Ti';
%targetGpuName = 'GeForce GTX 1070';

% You can get the list of GPU IDs
gpuidarray = getGpuIds(targetGpuName);
fprintf('Result of (%s)\n', targetGpuName);
fprintf('\t Number of devices: %d\n', length(gpuidarray));
disp(gpuidarray);

% If all the GPUs you have are the same and you are going to use them all,
% you don't add anything to the existing code.
% Your code will look like 
%    projections=Ax(head,geo,angles,'interpolated');
%    imgFDK=FDK(noise_projections,geo,angles);
%    imgOSSART=OS_SART(noise_projections,geo,angles,niter);
% or you can write like
%    gpuids = GpuIds();
%    projections=Ax(head,geo,angles,'interpolated', 'gpuids', gpuids);
%    imgFDK=FDK(noise_projections,geo,angles, 'gpuids', gpuids);
%    imgOSSART=OS_SART(noise_projections,geo,angles,niter, 'gpuids', gpuids);
%
% In the case where you have multiple types of GPUs installed or 
% you are going to limit the number of GPUs for a specific process,
% you can specify the GPUs using a GpuIds object.
% For example, if you have four 'GeForce RTX 2080 Ti' and a 'GeForce GTX 1060 6GB'
% and you are going to use all 'GeForce RTX 2080 Ti', use
%    gpuids = GpuIds('GeForce RTX 2080 Ti');
% If you use two of the GPUs, assuming whose device ID are 2 and 3, use
%    gpuids = GpuIds();
%    gpuids.devices = int32(2:3);

gpuids = GpuIds();
gpuids.devices = int32(1:10);
disp(gpuids);

gpuids2080 = GpuIds('GeForce RTX 2080 Ti');
disp(gpuids2080);


