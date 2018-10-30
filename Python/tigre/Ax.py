from _Ax import _Ax_ext
import numpy as np
import inspect
from matplotlib import pyplot as plt

def Ax(img, geo, angles,  krylov="interpolated"):
    #TODO: Does the raw image data match geo?
    if hasattr(geo,'check_geo'):
        geo.check_geo(angles)
    return _Ax_ext(img,geo,geo.angles,krylov,geo.mode)
