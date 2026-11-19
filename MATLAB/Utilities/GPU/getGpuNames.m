function names = getGpuNames()
    backend = tigreBackend();
    if backend == "cpu"
        names = {'CPU'};
        return
    elseif backend == "metal"
        if exist('tigre_metal_mex', 'file') == 3
            names = {tigre_metal_mex('deviceName')};
        else
            names = {'Apple Metal GPU'};
        end
        return
    end
    deviceCount = getGpuCount_mex();
    names = {''};
    for dev = 0:(deviceCount-1)
        name = getGpuName_mex(dev);
        names{1, dev+1}=name;
    end
end
