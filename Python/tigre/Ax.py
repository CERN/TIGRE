from _Ax import _Ax_ext

def Ax(img, geo, angles,  krylov="interpolated", mode="cone"):

    if hasattr(geo,'check_geo'):
        geo.check_geo(angles)
        print(geo.offOrigin.shape)
        print(geo.offOrigin)

    return _Ax_ext(img,geo,angles,krylov,mode)
