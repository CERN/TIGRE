from __future__ import division
import numpy as np
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
import time
import os
import math


class CGLS(IterativeReconAlg):
    __doc__ = (
        " CGLS_CBCT solves the CBCT problem using the conjugate gradient least\n"
        " squares\n"
        " \n"
        "  CGLS_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations.") + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.log_parameters = False
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

        if self.log_parameters:
            parameter_history = {}
            iterations = self.niter
            parameter_history['alpha'] = np.zeros(
                [iterations], dtype=np.float32)
            parameter_history['beta'] = np.zeros(
                [iterations], dtype=np.float32)
            parameter_history['gamma'] = np.zeros(
                [iterations], dtype=np.float32)
            parameter_history['q_norm'] = np.zeros(
                [iterations], dtype=np.float32)
            parameter_history['s_norm'] = np.zeros(
                [iterations], dtype=np.float32)
            self.parameter_history = parameter_history

        self.__r__ = self.proj - \
            Ax(self.res, self.geo, self.angles, 'ray-voxel')
        self.__p__ = Atb(self.__r__, self.geo, self.angles)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    def reinitialise_cgls(self):
        self.__r__ = self.proj - \
            Ax(self.res, self.geo, self.angles, 'ray-voxel')
        self.__p__ = Atb(self.__r__, self.geo, self.angles)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    # Overide
    def run_main_iter(self):
        self.l2l = np.zeros([self.niter], dtype=np.float32)
        avgtime = []
        for i in range(self.niter):
            if i == 0:
                if self.verbose:
                    print("CGLS Algorithm in progress.")
                toc = time.clock()
            if i == 1:
                tic = time.clock()
                if self.verbose:
                    print('Esitmated time until completetion (s): ' +
                      str((self.niter - 1) * (tic - toc)))
            avgtic = time.clock()
            q = tigre.Ax(self.__p__, self.geo, self.angles, 'ray-voxel')
            q_norm = np.linalg.norm(q)
            alpha = self.__gamma__ / (q_norm * q_norm)
            self.res += alpha * self.__p__
            avgtoc = time.clock()
            avgtime.append(abs(avgtic - avgtoc))
            for item in self.__dict__:
                if isinstance(getattr(self, item), np.ndarray):
                    if np.isnan(getattr(self, item)).any():
                        raise ValueError(
                            'nan found for ' + item + ' at iteraton ' + str(i))

            aux = self.proj - \
                tigre.Ax(self.res, self.geo, self.angles, 'ray-voxel')
            self.l2l[i] = np.linalg.norm(aux)
            if i > 0 and self.l2l[i] > self.l2l[i - 1]:
                if self.verbose:
                    print('re-initilization of CGLS called at iteration:' + str(i))
                if self.re_init_at_iteration + 1 == i:
                    if self.verbose:
                        print('Algorithm exited with two consecutive reinitializations.')
                    return self.res
                self.res -= alpha * self.__p__
                self.reinitialise_cgls()
                self.re_init_at_iteration = i

            self.__r__ -= alpha * q
            s = tigre.Atb(self.__r__, self.geo, self.angles)
            s_norm = np.linalg.norm(s)

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
        if self.verbose:
            print('Average time taken for each iteration for CGLS:' +
                  str(sum(avgtime) / len(avgtime)) + '(s)')


cgls = decorator(CGLS, name='cgls')
