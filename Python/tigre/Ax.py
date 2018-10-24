from _Ax import _Ax_ext
import numpy as np
import inspect
from matplotlib import pyplot as plt

def Ax(img, geo, angles,  krylov="interpolated", mode="cone"):

    if hasattr(geo,'check_geo'):
        geo.check_geo(angles,verbose=True)
    return _Ax_ext(img,geo,geo.angles,krylov,mode)
