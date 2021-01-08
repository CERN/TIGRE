
import _gpuUtils as gpuUtils

class GpuIds(object):
    """GpuIds
    To get all GPU names and Ids,
    ```
    gpuids = GpuIds()
    ```
    To get GPU Ids of specific name,
    ```
    gpuids = GpuIds('GeForce RTX 2080 Ti')
    ```
    """
    def __init__(self, nameGPU = None):
        if nameGPU is None:
            self.name = ""
        else:
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

# def getGpuCount():
#     # Returns number of installed GPUs.
#     return gpuUtils.getGpuCount()

def getGpuIds(gpuName):
    # Returns the GpuIds object, which contains the indexes of the device that whose name matches gpuName.
    gpuids = GpuIds(gpuName)
    return gpuids
