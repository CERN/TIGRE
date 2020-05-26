#%%
# 

isDebug = True
if isDebug:
    import sys
    sys.path.append('.')

from tigre.utilities import gpu

print("import _gpuUtils OK")


print("Querying installed GPU names")
listDeviceNames = gpu.getGpuNames()
print("\tDeviceCount: {}".format(len(listDeviceNames)))
print("\tNames      : {}".format(listDeviceNames))
print ("===================")


listQuery = [""]
listQuery.extend(listDeviceNames) # ["", 'GeForce GTX 1060 6GB', 'GeForce GTX 1080']

for nameGpu in listQuery:
    if nameGpu == "":
        print ("Querying All devices")
    else:
        print ("Querying {}".format(nameGpu))
    # gpuids = gpu.getGpuIdList(nameGpu)
    gpuids = gpu.getGpuIds(nameGpu)
    print("\tGpuIDs: {}".format(gpuids))
    print("\t # of IDs: {}".format(len(gpuids)))
    print("---")


