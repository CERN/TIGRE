from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator


class SART(IterativeReconAlg):
    __doc__ = ("solves the reconstruction problem\n"
               "using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "to the geometry described in GEO, using NITER iterations. \n")

    def __init__(self, proj, geo, angles, niter, **kwargs):
        kwargs.update(dict(blocksize=1))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


sart = decorator(SART, name='sart')


class SIRT(IterativeReconAlg):
    __doc__ = ("solves the reconstruction problem\n"
               "using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "to the geometry descrived in GEO, using NITER iterations.\n")

    def __init__(self, proj, geo, angles, niter, **kwargs):
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


sirt = decorator(SIRT, name='sirt')


class OS_SART(IterativeReconAlg):
    __doc__ = ( "solves the reconstruction problem\n"
                "using the projection data PROJ taken over ALPHA angles, corresponding\n"
                "to the geometry descrived in GEO, using NITER iterations.\n")

    def __init__(self, proj, geo, angles, niter, **kwargs):
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


ossart = decorator(OS_SART, name='ossart')
