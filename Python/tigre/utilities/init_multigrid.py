import numpy as np
import tigre
from numpy.lib.stride_tricks import as_strided
from scipy.interpolate import RegularGridInterpolator


def init_multigrid(proj, geo, alpha, alg):
    # WARNING: This takes a lot of memory!
    if alg == "SART":
        italg = tigre.algorithms.sart
    if alg == "SIRT":
        italg = tigre.algorithms.sirt
    finalsize = geo.nVoxel

    maxval = max(proj.ravel())
    minval = min(proj.ravel())

    # Start with 16 (raise this for larger images)
    geo.nVoxel = np.array([16, 16, 16])
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    if (geo.nVoxel > finalsize).all():
        return np.zeros(finalsize, dtype=np.float32)
    niter = 100
    initres = np.zeros(geo.nVoxel, dtype=np.float32)
    while (geo.nVoxel != finalsize).all():
        geo.dVoxel = geo.sVoxel / geo.nVoxel

        initres = italg(proj, geo, alpha, niter, init=initres, verbose=False)

        # get new dims(should be a way to do this more efficiently).

        geo.nVoxel = geo.nVoxel * 2
        geo.nVoxel[geo.nVoxel > finalsize] = finalsize[geo.nVoxel > finalsize]
        geo.dVoxel = geo.sVoxel / geo.nVoxel
        (x, y, z) = (
            np.linspace(minval, maxval, geo.nVoxel[0] / 2, dtype=np.float32),
            np.linspace(minval, maxval, geo.nVoxel[1] / 2, dtype=np.float32),
            np.linspace(minval, maxval, geo.nVoxel[2] / 2, dtype=np.float32),
        )

        # evaluate the function sart at the points xv,yv,zv

        xv, yv, zv = [
            tile_array(tile_array(x, 2), geo.nVoxel[0] ** 2),
            tile_array(tile_array(y, 2), geo.nVoxel[0] ** 2),
            tile_array(tile_array(x, 2), geo.nVoxel[0] ** 2),
        ]

        initres = RegularGridInterpolator((x, y, z), initres)(np.column_stack((xv, yv, zv)))
        initres = initres.reshape(geo.nVoxel)

    return initres


def tile_array(mat, b1):
    (r,) = mat.shape
    (rs,) = mat.strides
    x = as_strided(mat, (r, b1), (rs, 0))
    return x.reshape(r * b1)
