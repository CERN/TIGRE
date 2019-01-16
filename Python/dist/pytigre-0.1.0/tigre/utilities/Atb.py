from _Atb import _Atb_ext
import numpy as np
import copy
def Atb(projections, geo, angles,  krylov="matched"):

    if projections.dtype != np.float32:
        raise TypeError("Input data should be float32, not " + str(projections.dtype))
    if not np.isreal(projections).all():
        raise ValueError("Complex types not compatible for back projection.")
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)

    if projections.shape != (geox.angles.shape[0],geox.nDetector[0],geox.nDetector[1]):
        raise ValueError("Expected data shape: " + str((geox.angles.shape[0], geox.nDetector[0],geox.nDetector[1])) +
                         " not compatible with: " + str(projections.shape))

    return _Atb_ext(projections, geox, geox.angles, krylov,geox.mode)