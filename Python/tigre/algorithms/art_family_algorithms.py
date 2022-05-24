import copy

from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.im_3d_denoise import im3ddenoise



class SART(IterativeReconAlg):  
    __doc__ = (
        "SART solves Cone Beam CT image reconstruction using \n"
        "Simultaneous Algebraic Reconstruction Technique algorithm\n"
        "SART(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry described in GEO, using NITER iterations. \n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to 1, please use an OS version of the algorithm for blocksize > 1')
        kwargs.update(dict(blocksize=1))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


sart = decorator(SART, name="sart")


class SIRT(IterativeReconAlg):  
    __doc__ = (
        "SIRT solves Cone Beam CT image reconstruction using \n"
        "Simultaneous Iterrative Reconxtructive Technique algorithm\n"
        "SIRT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry descrived in GEO, using NITER iterations.\n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to {}, please do not specify blocksize for this algorithm'.format(angles.shape[0]))
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


sirt = decorator(SIRT, name="sirt")


class OS_SART(IterativeReconAlg):  
    __doc__ = (
        "OS_SART solves Cone Beam CT image reconstruction using Oriented Subsets\n"
        "Simultaneous Algebraic Reconxtruction Techique algorithm\n"
        "OS_SART(PROJ,GEO,ALPHA,NITER,BLOCKSIZE=20) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry descrived in GEO, using NITER iterations.\n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        
        self.blocksize = 20 if 'blocksize' not in kwargs else kwargs["blocksize"]       
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


ossart = decorator(OS_SART, name="ossart")


class SART_TV(IterativeReconAlg):  
    __doc__ = (
        "SART_TV solves Cone Beam CT image reconstruction using Simultaneous \n"
        "Algebraic Reconstruction Technique with TV regularization algorithm\n"
        "SART_TV(PROJ,GEO,ALPHA,NITER,TVLAMBDA=50,TVITER=50) solves the reconstruction\n"
        "problem using the projection data PROJ taken over ALPHA angles\n"
        "corresponding to the geometry described in GEO, using NITER iterations. \n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        
        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to 1, please use an OS version of the algorithm for blocksize > 1')
        kwargs.update(dict(blocksize=1))
        self.tvlambda = 50 if 'tvlambda' not in kwargs else kwargs['tvlambda']
        self.tviter = 50 if 'tviter' not in kwargs else kwargs['tviter']
        # these two settings work well for nVoxel=[254,254,254]

        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    # Override
    def run_main_iter(self):
        """
        Goes through the main iteration for the given configuration.
        :return: None
        """
        Quameasopts = self.Quameasopts

        for i in range(self.niter):

            res_prev = None
            if Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            if self.verbose:
                self._estimate_time_until_completion(i)

            getattr(self, self.dataminimizing)()
            # print("run_main_iter: gpuids = {}", self.gpuids)
            self.res = im3ddenoise(self.res, self.tviter, self.tvlambda, self.gpuids)
            if Quameasopts is not None:
                self.error_measurement(res_prev, i)


sart_tv = decorator(SART_TV, name="sart_tv")


class OSSART_TV(IterativeReconAlg):  
    __doc__ = (
        "OSSART_TV solves Cone Beam CT image reconstruction using Oriented Subsets\n"
        "Simultaneous Algebraic Reconxtruction Technique with TV regularization algorithm\n"
        "OSSART_TV(PROJ,GEO,ALPHA,NITER,BLOCKSIZE=20,TVLAMBDA=50,TVITER=50) \n"
        "solves the reconstruction problem using the projection data PROJ taken\n"
        "over ALPHA angles, corresponding to the geometry descrived in GEO,\n"
        "using NITER iterations.\n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        
        self.blocksize = 20 if 'blocksize' not in kwargs else kwargs['blocksize']
        self.tvlambda = 50 if 'tvlambda' not in kwargs else kwargs['tvlambda']
        self.tviter = 50 if 'tviter' not in kwargs else kwargs['tviter']
        # these two settings work well for nVoxel=[254,254,254]

        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


ossart_tv = decorator(OSSART_TV, name="ossart_tv")
