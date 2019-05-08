from _Ax import _Ax_ext
import numpy as np
import copy

def Ax(img, geo, angles,  krylov="ray-voxel"):

    if img.dtype != np.float32:
        raise TypeError("Input data should be float32, not "+ str(img.dtype))
    if not np.isreal(img).all():
        raise ValueError("Complex types not compatible for projection.")
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    """
    Here we cast all values in geo to single point precision float. This way we 
    know what behaviour to expect from pytigre to Cuda and can change 
    single parameters accordingly. 
    """
    geox.cast_to_single()
    geox.checknans()

    if all(img.shape != geox.nVoxel):
        raise ValueError("Input data should be of shape geo.nVoxel: "+ str(geox.nVoxel) +
                         " not:" + str(img.shape))

    return _Ax_ext(img, geox, geox.angles, krylov, geox.mode)
