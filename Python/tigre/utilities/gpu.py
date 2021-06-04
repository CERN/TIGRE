"""
1. Check the names of installed GPUs.
   names = getGpuNames()
   print(names)
2. Select the name of the GPUs that you want your program run on.
   gpuids = getGpuIds('GeForce RTX 2080 Ti')
3. Pass the gpuids instance to Tigre functions.
"""
import _gpuUtils as gpuUtils


class GpuIds(object):
    """
    A class that holds the IDs and their (common) name
    ```
    gpuids = GpuIds('GeForce RTX 2080 Ti')
    ```
    If name is ""/None, IDs of all GPUs are hold.
    ```
    gpuids = GpuIds()
    ```
    """

    def __init__(self, nameGPU=None):  # noqa: N803
        if nameGPU is None:
            self.name = ""
        else:
            self.name = nameGPU
        self.devices = gpuUtils.getGpuIdList(self.name)

    def __len__(self):
        return len(self.devices)

    def __str__(self):
        dictTemp = {
            "name": self.name,
            "devices": self.devices,
        }
        return dictTemp.__str__()


def getGpuNames():
    """Returns a list of all installed GPUs."""
    return gpuUtils.getGpuNames()


# def getGpuCount():
#     """Returns the number of installed GPUs."""
#     return gpuUtils.getGpuCount()


def getGpuIds(gpuName=None):  # noqa: N803
    """Returns the GpuIds object, which contains the IDs of devices whose name matches gpuName."""
    gpuids = GpuIds(gpuName)
    return gpuids
