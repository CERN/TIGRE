import copy

import numpy as np
from _Ax import _Ax_ext

from .gpu import GpuIds


def Ax(img, geo, angles, projection_type="Siddon", **kwargs):
    if img.dtype != np.float32:
        raise TypeError("Input data should be float32, not " + str(img.dtype))
    if not np.isreal(img).all():
        raise ValueError("Complex types not compatible for projection.")
    if any(img.shape != geo.nVoxel):
        raise ValueError(
            "Input data should be of shape geo.nVoxel: "
            + str(geo.nVoxel)
            + " not:"
            + str(img.shape)
        )
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    """
    Here we cast all values in geo to single point precision float. This way we know what behavior
    to expect from pytigre to Cuda and can change single parameters accordingly.
    """
    geox.cast_to_single()
    # geox.checknans()

    if kwargs.get("gpuids", None) is None:
        gpuids = GpuIds()
    else:
        gpuids = kwargs["gpuids"]
    # if we have more GPUs than projections to compute, reduce the amount of GPUs. 
    if len(angles) < len(gpuids):
        gpuids.devices = list(gpuids.devices[0:len(angles)])

    return _Ax_ext(img, geox, geox.angles, projection_type, geox.mode, gpuids=gpuids)
