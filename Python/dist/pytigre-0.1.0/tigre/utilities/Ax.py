from _Ax import _Ax_ext
import numpy as np
import copy

def Ax(img, geo, angles,  krylov="interpolated"):

    if img.dtype != np.float32:
        raise TypeError("Input data should be float32, not "+ str(img.dtype))
    if not np.isreal(img).all():
        raise ValueError("Complex types not compatible for projection.")
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)


    if all(img.shape != geox.nVoxel):
        raise ValueError("Input data should be of shape geo.nVoxel: "+ str(geox.nVoxel) +
                         " not:" + str(img.shape))

    return _Ax_ext(img,geox,geox.angles,krylov,geox.mode)
