from _Atb import _Atb_ext
import copy
def Atb(projections, geo, angles,  krylov="matched"):

    geo.check_geo(angles)

    return _Atb_ext(projections, geo, geo.angles, krylov,geo.mode)