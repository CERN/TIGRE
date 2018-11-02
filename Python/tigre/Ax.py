from _Ax import _Ax_ext
import copy

def Ax(img, geo, angles,  krylov="interpolated"):
    #TODO: Does the raw image data match geo?

    geo.check_geo(angles)

    return _Ax_ext(img,geo,geo.angles,krylov,geo.mode)
