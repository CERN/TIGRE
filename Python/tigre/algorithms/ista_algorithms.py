from __future__ import division
import numpy as np
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
import time
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

        # Dont precompute W and V
        kwargs.update(dict(W=None,
                           V=None,
                           ))
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        self.lmbda = 0.1

        self.bigl = 2.e8
        self.bm = 1./self.bigl

    # overide update_image from iterative recon alg to remove W.
    def update_image(self, geo, angle, iteration):
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
        self.res += self.bm * 2 * tigre.Atb((self.proj[self.angle_index[iteration]]
                                                                     - tigre.Ax(self.res, geo, angle, 'interpolated')),
                                                                    geo, angle, 'matched')
    def run_main_iter(self):
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
                        print(str(self.name).upper() + ' ' + "algorithm in progress.")
                        toc = time.clock()
                    if i == 1:
                        tic = time.clock()
                        print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
                getattr(self, self.dataminimizing)()

                self.error_measurement(res_prev, i)


fista = decorator(FISTA)
