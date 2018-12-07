from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator


class SART(IterativeReconAlg):
    __doc__ = ("SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
               "Simultaneous Algebraic Reconstruction Techique algorithm\n"
               "SART(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem\n"
               "using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "to the geometry described in GEO, using NITER iterations. \n") + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        kwargs.update(dict(blocksize=1))
        IterativeReconAlg.__init__(self,proj, geo, angles, niter, **kwargs)



sart = decorator(SART, name='sart')


class SIRT(IterativeReconAlg):
    __doc__ = ("SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
               "Simultaneous Algebraic Reconxtruction Techique algorithm\n"

               "SIRT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem\n"
               "using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "to the geometry descrived in GEO, using NITER iterations.\n") + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self,proj, geo, angles, niter, **kwargs)



sirt = decorator(SIRT, name='sirt')


class OS_SART(IterativeReconAlg):
    __doc__ = ("SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
               "Simultaneous Algebraic Reconxtruction Techique algorithm\n"

               "OS_SART(PROJ,GEO,ALPHA,NITER,BLOCKSIZE=20) solves the reconstruction problem\n"
               "using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "to the geometry descrived in GEO, using NITER iterations.\n") + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        IterativeReconAlg.__init__(self,proj, geo, angles, niter, **kwargs)


ossart = decorator(OS_SART, name='ossart')
