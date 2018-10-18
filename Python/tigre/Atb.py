from _Atb import _Atb_ext

def Atb(projections, geo, angles,  krylov="matched", mode="cone"):
    if hasattr(geo,'check_geo'):
        geo.check_geo(angles)
        return _Atb_ext(projections, geo, angles, krylov,mode)