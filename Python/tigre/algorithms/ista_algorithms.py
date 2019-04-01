from __future__ import division
import numpy as np
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
import time
from tigre.utilities.im_3d_denoise import im3ddenoise
from tigre.algorithms.single_pass_algorithms import FDK
import copy


class FISTA(IterativeReconAlg):
    __doc__ = (
                  """
                  
                  :param specific kwargs:
      
                      :keyword largest_eigen_value: Estimated largest eigenvalue 
                      squared of the matrix A in the equations Ax or At - b
                  """) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        # Dont precompute W and V
        kwargs.update(dict(W=None,
                           V=None,
                           ))
        kwargs.update(dict(blocksize=angles.shape[0]))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        self.lmbda = 0.1
        if 'largest_singular_value' not in kwargs:
            self.__L__ = 2.e4
        else:
            self.__L__ = kwargs['largest_singular_value']
        self.__t__ = 1
        self.__bm__ = 1. / self.__L__

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
        self.res += self.__bm__ * 2 * tigre.Atb((self.proj[self.angle_index[iteration]]
                                                 - tigre.Ax(self.res, geo, angle, 'interpolated')),
                                                geo, angle, 'matched')

    def run_main_iter(self):
        """
        Goes through the main iteration for the given configuration.
        :return: None
        """
        t = self.__t__
        Quameasopts = self.Quameasopts
        x_rec = copy.deepcopy(self.res)
        lambdaForTv = 2 * self.__bm__ * self.lmbda
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

            x_rec_old = copy.deepcopy(x_rec)
            x_rec = im3ddenoise(self.res, 20, 1. / lambdaForTv)
            t_old = t
            t = (1 + np.sqrt(1 + 4 * t ** 2)) / 2
            self.res = x_rec + (t_old - 1) / t * (x_rec - x_rec_old)

            self.error_measurement(res_prev, i)


fista = decorator(FISTA,name='FISTA')


class ISTA(FISTA):
    __doc__ = FISTA.__doc__

    def __int__(self, proj, geo, angles, niter, **kwargs):
        FISTA.__init__(self, proj, geo, angles, niter, **kwargs)

    def run_main_iter(self):

        """
        Goes through the main iteration for the given configuration.
        :return: None
        """
        Quameasopts = self.Quameasopts
        lambdaForTv = 2 * self.__bm__ * self.lmbda
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

            self.res = im3ddenoise(self.res, 20, 1. / lambdaForTv)

            self.error_measurement(res_prev, i)

ista = decorator(ISTA,name='ISTA')