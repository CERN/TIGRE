from __future__ import division

import copy
import numpy as np
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Atb import Atb
from tigre.utilities.Ax import Ax


class MLEM(IterativeReconAlg):  # noqa: D101
    __doc__ = (
        " MLEM solves the CBCT problem using the maximum likelihood expectation maximization\n"
        " algorithm\n"
        " \n"
        "  MLEM(PROJ,GEO,ANGLES,NITER,INIT) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations. INIT specifies\n"
        "  starting image, default: None (flat image value=1)"
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

        if self.init is None:
            self.res += 1.0

        self.W = Atb(np.ones(proj.shape, dtype=np.float32), geo, angles, backprojection_type="matched", gpuids=self.gpuids)
        self.W[self.W <= 0.0] = np.inf

    # override
    def run_main_iter(self):
        self.res[self.res < 0.0] = 0.0
        for i in range(self.niter):
            res_prev = copy.deepcopy(self.res) if self.Quameasopts is not None else None
            if self.verbose:
                self._estimate_time_until_completion(i)

            den = Ax(self.res, self.geo, self.angles, "interpolated", gpuids=self.gpuids)
            den[den == 0.0] = np.inf
            auxmlem = self.proj / den

            # update
            img = Atb(auxmlem, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids) / self.W  

            self.res = self.res * img
            self.res[self.res < 0.0] = 0.0
            # with computel2=True, records the KL divergence (see error_measurement)
            self.error_measurement(res_prev, i)

    # override: MLEM's data error is the Kullback-Leibler divergence, not the L2 residual
    def error_measurement(self, res_prev, iter):
        computel2 = self.computel2
        self.computel2 = False  # the base class records the quality measures only
        try:
            IterativeReconAlg.error_measurement(self, res_prev, iter)
        finally:
            self.computel2 = computel2
        if computel2 and iter < self.l2l.shape[1]:
            est = Ax(self.res, self.geo, self.angles, "interpolated", gpuids=self.gpuids)
            self.l2l[0, iter] = kl_divergence(self.proj, est)


def kl_divergence(proj, est):
    """Kullback-Leibler divergence of the estimate ``est = Ax`` from the data ``proj = b``.

    ``sum(b * log(b / Ax) - b + Ax)``: the Poisson negative log-likelihood up to
    a constant, which MLEM decreases. Pixels with ``b <= 0`` count as ``Ax``
    (``0 * log(0) = 0``). Rays with ``Ax = 0`` are left out: the multiplicative
    update cannot change them.
    """
    est = np.asarray(est)
    proj = np.asarray(proj)
    use = est > 0
    a = est[use].astype(np.float64)
    b = proj[use].astype(np.float64)
    pos = b > 0
    return float(a.sum() - b[pos].sum() + np.sum(b[pos] * np.log(b[pos] / a[pos])))


mlem = decorator(MLEM, name="mlem")
