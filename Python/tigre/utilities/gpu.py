
import _gpuUtils as gpuUtils

class GpuIds(object):
    def __init__(self):
        self.name = ""
        self.devices = []
    def setup(self, nameGPU):
        self.name = nameGPU
        self.devices = gpuUtils.getGpuIdList(nameGPU)
    def __len__(self):
        return len(self.devices)
    def __str__(self):
        dictTemp = {
            "name": self.name,
            "devices": self.devices,
        }
        return dictTemp.__str__()


def getGpuNames():
    # Returns a list of all installed GPUs.
    return gpuUtils.getGpuNames()

def getGpuIdList(gpuName):
    # Returns the list of index of the device that whose name maches gpuName.
    return gpuUtils.getGpuIdList()

def getGpuIds(gpuName):
    # Returns the GpuIds object, which contains the indexes of the device that whose name maches gpuName.
    gpuids = GpuIds()
    gpuids.setup(gpuName)
    return gpuids
