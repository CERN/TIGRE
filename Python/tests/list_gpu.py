#%%
isDebug = True
if isDebug:
    import os
    import sys
    dirsup = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    sys.path.insert(0, dirsup)

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
    gpuids = gpu.getGpuIds(nameGpu)
    print("\tGpuIDs: {}".format(gpuids))
    print("\t # of IDs: {}".format(len(gpuids)))
    print("---")




# %%
