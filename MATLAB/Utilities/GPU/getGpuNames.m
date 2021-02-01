function names = getGpuNames()
    deviceCount = getGpuCount_mex();
    names = {''};
    for dev = 0:(deviceCount-1)
        name = getGpuName_mex(dev);
        names{1, dev+1}=name;
    end
end
