import tigre
from tigre.utilities import gpu
# List the names of installed GPUs
print("Querying installed GPU names")
listDeviceNames = gpu.getGpuNames()  # noqa: N816
print("\tDeviceCount: {}".format(len(listDeviceNames)))
print("\tNames      : {}".format(listDeviceNames))
print("===================")

# Choose one of them
# targetGpuName=''  # noqa: N816
targetGpuName = "GeForce GTX 1060 6GB"  # noqa: N816
# targetGpuName = 'GeForce RTX 2080 Ti'  # noqa: N816
# targetGpuName = 'GeForce GTX 1070'  # noqa: N816

# You can get the list of GPU IDs
gpuids = gpu.getGpuIds(targetGpuName)
print("Result of ({})".format(targetGpuName))
print("\t Number of devices: {}".format(len(gpuids)))
print("\t gpuids: {}".format(gpuids))


# If all the GPUs you have are the same and you are going to use them all, you don't add anything to
# the existing code.
# Your code will look like:

# projections=tigre.Ax(head,geo,angles)
# imgFDK=tigre.fdk(noise_projections,geo,angles)
# imgOSSART=tigre.fdk.os_sart(noise_projections,geo,angles,niter)

# or you can write like:

# projections=tigre.Ax(head,geo,angles,gpuids=gpuids)
# imgFDK=tigre.fdk(noise_projections,geo,angles,gpuids=gpuids)
# imgOSSART=tigre.fdk.os_sart(noise_projections,geo,angles,niter,gpuids=gpuids)


# In the case where you have multiple types of GPUs installed or you are going to limit the number
# of GPUs for a specific process, you can specify the GPUs using a GpuIds object.
# For example, if you have four 'GeForce RTX 2080 Ti' and a 'GeForce GTX 1060 6GB' and you are going
# to use all 'GeForce RTX 2080 Ti', use
#    gpuids = GpuIds('GeForce RTX 2080 Ti');
# If you use two of the GPUs, assuming whose device ID are 2 and 3, use
#    gpuids = GpuIds();
#    gpuids.devices = int32(2:3);

gpuids = gpu.GpuIds()
gpuids.devices = list(range(0, 10))
print(gpuids)

gpuids2080 = gpu.GpuIds("GeForce RTX 2080 Ti")
print(gpuids2080)


