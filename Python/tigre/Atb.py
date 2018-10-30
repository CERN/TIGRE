from _Atb import _Atb_ext
import numpy as np
import inspect

def Atb(projections, geo, angles,  krylov="matched"):
    if hasattr(geo,'check_geo'):
        geo.check_geo(angles)
        return _Atb_ext(projections, geo, geo.angles, krylov,geo.mode)