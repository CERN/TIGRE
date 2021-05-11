import copy
import time

from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.im_3d_denoise import im3ddenoise


if hasattr(time, "perf_counter"):
    default_timer = time.perf_counter
else:
    default_timer = time.clock


class SART(IterativeReconAlg):  # noqa: D101
    __doc__ = (
        "SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
        "Simultaneous Algebraic Reconstruction Techique algorithm\n"
        "SART(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry described in GEO, using NITER iterations. \n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        kwargs.update(dict(blocksize=1))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


sart = decorator(SART, name="sart")


class SIRT(IterativeReconAlg):  # noqa: D101
    __doc__ = (
        "SIRT_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
        "Simultaneous Algebraic Reconxtruction Techique algorithm\n"
        "SIRT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry descrived in GEO, using NITER iterations.\n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


sirt = decorator(SIRT, name="sirt")


class OS_SART(IterativeReconAlg):  # noqa: D101, N801
    __doc__ = (
        "OS_SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
        "Simultaneous Algebraic Reconxtruction Techique algorithm\n"
        "OS_SART(PROJ,GEO,ALPHA,NITER,BLOCKSIZE=20) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry descrived in GEO, using NITER iterations.\n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)


ossart = decorator(OS_SART, name="ossart")


class OS_SART_TV(IterativeReconAlg):  # noqa: D101, N801
    __doc__ = (
        "OS_SART_TV_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets\n"
        "Simultaneous Algebraic Reconxtruction Techique algorithm\n"
        "OS_SART_TV(PROJ,GEO,ALPHA,NITER,BLOCKSIZE=20) solves the reconstruction problem\n"
        "using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "to the geometry descrived in GEO, using NITER iterations.\n"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        if "tvlambda" not in kwargs:
            kwargs.update(dict(tvlambda=50))
        if "tviter" not in kwargs:
            kwargs.update(dict(tviter=50))
        # these two settings work well for nVoxel=[254,254,254]

        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

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
                if i == 0:
                    print(str(self.name).upper() + " " + "algorithm in progress.")
                    toc = default_timer()
                if i == 1:
                    tic = default_timer()

                    remaining_time = (self.niter - 1) * (tic - toc)
                    seconds = int(remaining_time)
                    print(
                        "Estimated time until completion : "
                        + time.strftime("%H:%M:%S", time.gmtime(seconds))
                    )

            getattr(self, self.dataminimizing)()
            # print("run_main_iter: gpuids = {}", self.gpuids)
            self.res = im3ddenoise(self.res, self.tviter, self.tvlambda, self.gpuids)
            self.error_measurement(res_prev, i)


ossart_tv = decorator(OS_SART_TV, name="ossart_tv")
