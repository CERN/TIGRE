function gpuids = getGpuIds(gpuname)
    deviceCount = getGpuCount_mex();
    gpuids = int32(0:(deviceCount-1));
    for idx = 1:deviceCount
        name = getGpuName_mex(gpuids(idx));
        if ~strcmp(name, gpuname)
            gpuids(idx) = -1;
        end
    end
    gpuids = gpuids(gpuids>=0);
end
