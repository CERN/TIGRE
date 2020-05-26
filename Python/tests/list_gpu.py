#%%
# 

# import _Ax
# print("_AX ok")
# import _Atb
# print("_Atb ok")

import _gpuUtils as gpu
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
    listDevices = gpu.getGpuIdList(nameGpu)
    print("\tlistDevices: {}".format(listDevices))
    print("---")



# %%
