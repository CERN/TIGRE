"""
1. Check the names of installed GPUs.
   names = getGpuNames()
   print(names)
2. Select the name of the GPUs that you want your program run on.
   gpuids = getGpuIds('GeForce RTX 2080 Ti')
3. Pass the gpuids instance to Tigre functions.
"""
from tigre.utilities.cuda_interface import _gpuUtils as gpuUtils
import copy
import importlib

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

    @classmethod
    def _from_devices(cls, devices):
        """GpuIds for these device numbers, without querying the GPUs."""
        ids = cls.__new__(cls)
        ids.name = ""
        ids.devices = [int(d) for d in devices]
        return ids

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


class _NotBuilt(object):
    # Stands in for an extension that the SYCL build does not have yet.
    def __init__(self, module, error):
        self._module = module
        self._error = error

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        module, error = self._module, self._error

        def missing(*args, **kwargs):
            raise ImportError(
                "{} is not in the SYCL build of TIGRE yet ({})".format(module, error), name=module
            )

        return missing


def _import_extension(module):
    """Import one of TIGRE's compiled extensions.

    The SYCL build does not have every extension yet (SYCL_PENDING in setup.py).
    There, a missing one is replaced by a stand-in whose functions raise the
    ImportError when they are called, so the rest of TIGRE can still be imported.
    """
    name = "tigre.utilities.cuda_interface." + module
    try:
        return importlib.import_module(name)
    except ImportError as err:
        if getattr(gpuUtils, "getBackend", lambda: "cuda")() != "sycl":
            raise
        return _NotBuilt(name, err)
