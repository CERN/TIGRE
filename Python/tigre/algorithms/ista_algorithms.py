from __future__ import division
import numpy as np
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.im3Dnorm import im3DNORM
from tigre.algorithms.single_pass_algorithms import FDK
import copy

class FISTA(IterativeReconAlg):
    def __init__(self, proj, geo, angles, niter, **kwargs):

        """

        :param proj:
        :param geo:
        :param angles:
        :param niter:
        :param kwargs:
        """

        # Dont precompute W, set V = ones.
        kwargs.update(dict(W=None,
                           V=np.ones((angles.shape[0]), dtype=np.float32)
                           ))
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        self.lmbda = 0.1
        self.bigl = 2.e5
        self.bm = 1./self.bigl

    # overide gradient descent from iterative recon alg to remove W.
    def gradient_descent(self, geo, angle, iteration):
        """
        VERBOSE:
         for j in range(angleblocks):
             angle = np.array([alpha[j]], dtype=np.float32)
             proj_err = proj[angle_index[j]] - Ax(res, geo, angle, 'ray-voxel')
             backprj = Atb(proj_err, geo, angle, 'FDK')
             res += backprj
             res[res<0]=0

        :return: None
        """
        self.res += self.lmbda * 1. / self.V[iteration] * tigre.Atb((self.proj[self.angle_index[iteration]]
                                                                     - tigre.Ax(self.res, geo, angle, 'interpolated')),
                                                                    geo, angle, 'FDK')
    def run_main_iter(self):
        for i in range(self.niter):
            y_rec = copy.deepcopy(self.res)
            self.art_data_minimizing()
            y_rec += self.bm*2*self.res


fista = decorator(FISTA)
