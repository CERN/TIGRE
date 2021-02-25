from __future__ import division
import numpy as np
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
import time


if hasattr(time, "perf_counter"):
    default_timer = time.perf_counter
else:
    default_timer = time.clock


class MLEM(IterativeReconAlg):
    __doc__ = (
        " MLEM_CBCT solves the CBCT problem using the maximum likelihood expectation maximization\n"
        " algorithm\n"
        " \n"
        "  MLEM_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

        if self.init is None:
            self.res += 1.0

        self.W = Atb(np.ones(proj.shape, dtype=np.float32), geo, angles, gpuids=self.gpuids)
        self.W[self.W <= 0.0] = np.inf

    # Overide
    def run_main_iter(self):
        self.res[self.res < 0.0] = 0.0
        for i in range(self.niter):
            if i == 0:
                if self.verbose:
                    print("MLEM Algorithm in progress.")
                toc = default_timer()
            if i == 1:
                tic = default_timer()
                if self.verbose:
                    print(
                        "Esitmated time until completetion (s): {:.4f}".format(
                            (self.niter - 1) * (tic - toc)
                        )
                    )

            #            tic = time.process_time()
            den = Ax(self.res, self.geo, self.angles, "interpolated", gpuids=self.gpuids)
            # toc = time.process_time()
            # print('Ax time: {}'.format(toc-tic))
            den[den == 0.0] = np.inf
            auxmlem = self.proj / den
            # auxmlem[auxmlem == np.nan] = 0.
            # auxmlem[auxmlem == np.inf] = 0.

            # update
            # tic = time.process_time()
            img = Atb(auxmlem, self.geo, self.angles, gpuids=self.gpuids) / self.W
            # toc = time.process_time()
            # print('Atb time: {}'.format(toc-tic))
            # img[img == np.nan] = 0.
            # img[img == np.inf] = 0.

            self.res = self.res * img
            self.res[self.res < 0.0] = 0.0


mlem = decorator(MLEM, name="mlem")
