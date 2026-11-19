function gpuids = getGpuIds(gpuname)
    if tigreBackend() ~= "cuda"
        names = getGpuNames();
        if any(strcmp(names, gpuname))
            gpuids = int32(0);
        else
            gpuids = int32.empty(1, 0);
        end
        return
    end
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
