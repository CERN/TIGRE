import numpy as np
from _tv_proximal import tvdenoise

from .gpu import GpuIds


def im3ddenoise(img, iter=50, lmbda=15.0, gpuids=None):
    if img.ndim != 3:
        raise ValueError(
            f"Expected a 3D array (nZ, nY, nX), got shape {img.shape}. "
            f"For 2D denoising, use tigre.utilities.im2ddenoise()."
        )
    imgmin = np.amin(img.ravel())
    img = img - imgmin
    imgmax = np.amax(img.ravel())
    img = img / imgmax

    if gpuids is None:
        gpuids = GpuIds()

    img = tvdenoise(img, iter, lmbda, gpuids)

    img *= imgmax
    img += imgmin

    return img
