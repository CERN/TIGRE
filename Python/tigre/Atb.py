from _Atb import _Atb_ext
import numpy as np
import inspect

def Atb(projections, geo, angles,  krylov="matched", mode="cone"):
    if hasattr(geo,'check_geo'):
        geo.check_geo(angles)
        print(geo.mode + ': angles: ' + str(geo.angles.shape))
        return _Atb_ext(projections, geo, geo.angles, krylov,mode)