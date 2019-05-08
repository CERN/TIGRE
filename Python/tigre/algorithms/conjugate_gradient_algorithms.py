from __future__ import division
import numpy as np
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.Ax import Ax
from tigre.Atb import Atb
import time
import os
import math


class CGLS(IterativeReconAlg):
    __doc__ = (" CGLS_CBCT solves the CBCT problem using the conjugate gradient least\n"
               " squares\n"
               " \n"
               "  CGLS_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
               "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "  to the geometry descrived in GEO, using NITER iterations.")

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.log_parameters = False
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

        if self.log_parameters:
            parameter_history = {}
            iterations = self.niter
            parameter_history['alpha'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['beta'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['gamma'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['q_norm'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['s_norm'] = np.zeros([iterations], dtype=np.float32)
            self.parameter_history = parameter_history

        self.__r__ = self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel')
        self.__p__ = Atb(self.__r__, self.geo, self.angles)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    def reinitialise_cgls(self):
        self.__r__ = self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel')
        self.__p__ = Atb(self.__r__, self.geo, self.angles)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    # Overide
    def run_main_iter(self):
        self.l2l = np.zeros([self.niter], dtype=np.float32)
        for i in range(self.niter):
            if i == 0:
                print("CGLS Algorithm in progress.")
                toc = time.clock()
            if i == 1:
                tic = time.clock()
                print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
            q = Ax(self.__p__, self.geo, self.angles, 'ray-voxel')
            q_norm = np.linalg.norm(q.ravel(), 2)
            alpha = self.__gamma__ / (q_norm * q_norm)
            self.res += alpha * self.__p__

            for item in self.__dict__:
                if type(getattr(self, item)) == np.ndarray:
                    if np.isnan(getattr(self, item)).any():
                        raise ValueError('nan found for ' + item + ' at iteraton ' + str(i))

            aux = self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel')
            self.l2l[i] = np.linalg.norm(aux.ravel(), 2)
            if i > 0 and self.l2l[i] > self.l2l[i - 1]:
                print('re-initialization was called at iter:' + str(i))
                self.res -= alpha * self.__p__
                self.reinitialise_cgls()

            self.__r__ -= alpha * q
            s = Atb(self.__r__, self.geo, self.angles)
            s_norm = np.linalg.norm(s.ravel(), 2)

            gamma1 = s_norm * s_norm
            beta = gamma1 / self.__gamma__
            if self.log_parameters:
                self.parameter_history['alpha'][i] = alpha
                self.parameter_history['beta'][i] = beta
                self.parameter_history['gamma'][i] = self.__gamma__
                self.parameter_history['q_norm'][i] = q_norm
                self.parameter_history['s_norm'][i] = s_norm

            self.__gamma__ = gamma1
            self.__p__ = s + beta * self.__p__


cgls = decorator(CGLS, name='cgls')
