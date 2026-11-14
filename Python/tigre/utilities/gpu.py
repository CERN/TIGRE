"""
1. Check the names of installed GPUs.
   names = getGpuNames()
   print(names)
2. Select the name of the GPUs that you want your program run on.
   gpuids = getGpuIds('GeForce RTX 2080 Ti')
3. Pass the gpuids instance to Tigre functions.
"""
import _gpuUtils as gpuUtils
import copy

class GpuIds(object):
    """A set of GPUs for TIGRE to use.

    Pass it to any operator or algorithm with ``gpuids=``.

    Parameters
    ----------
    nameGPU : str, optional
        Only GPUs with exactly this name (see :func:`getGpuNames`). Default:
        all GPUs.

    Attributes
    ----------
    devices : list of int
        CUDA device numbers. Set it to choose devices by number.
    name : str
        The name the GPUs were selected by.

    Examples
    --------
    >>> from tigre.utilities.gpu import GpuIds, getGpuNames
    >>> getGpuNames()
    ['NVIDIA RTX A6000', 'NVIDIA RTX A6000']
    >>> gpuids = GpuIds()          # all GPUs
    >>> gpuids.devices = [0]       # only the first one
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
    def __getitem__(self, key):
        single_gpu = copy.deepcopy(self)
        single_gpu.devices = [self.devices[key]]
        return single_gpu

def getGpuNames():
    """Names of the installed GPUs, one per device.

    Returns
    -------
    list of str
    """
    return gpuUtils.getGpuNames()


# def getGpuCount():
#     """Returns the number of installed GPUs."""
#     return gpuUtils.getGpuCount()


def getGpuIds(gpuName=None):  # noqa: N803
    """GPUs whose name matches ``gpuName``; the same as ``GpuIds(gpuName)``.

    Returns
    -------
    GpuIds
    """
    gpuids = GpuIds(gpuName)
    return gpuids
