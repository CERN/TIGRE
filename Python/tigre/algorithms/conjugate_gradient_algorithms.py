from __future__ import division
import numpy as np
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
import time
class CGLS(IterativeReconAlg):
    __doc__ = (" CGLS_CBCT solves the CBCT problem using the conjugate gradient least\n"
               " squares\n"
               " \n"
               "  CGLS_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
               "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
               "  to the geometry descrived in GEO, using NITER iterations.") +IterativeReconAlg.__doc__

    def __init__(self,proj,geo,angles,niter,**kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None,V=None))
        self.log_parameters = False
        # Avoid typo checking
        IterativeReconAlg.__init__(self,proj,geo,angles,niter,**kwargs)

        self.initialise_cgls()


        if self.log_parameters:
            parameter_history = {}
            iterations = self.niter
            parameter_history['alpha'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['beta'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['gamma'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['q_norm'] = np.zeros([iterations], dtype=np.float32)
            parameter_history['s_norm'] = np.zeros([iterations], dtype=np.float32)
            self.parameter_history = parameter_history

    def initialise_cgls(self):
        self._r = self.proj - Ax(self.res,self.geo,self.angles,'ray-voxel')
        self._p = Atb(self._r,self.geo,self.angles)
        p_norm = np.linalg.norm(self._p.ravel(), 2)
        self._gamma = p_norm * p_norm

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
            q = Ax(self._p, self.geo, self.angles, 'ray-voxel')
            q_norm = np.linalg.norm(q.ravel(), 2)
            alpha = self._gamma / (q_norm * q_norm)
            self.res += alpha * self._p

            aux = self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel')
            self.l2l[i] = np.linalg.norm(aux.ravel(), 2)

            if i > 0 and self.l2l[i] > self.l2l[i - 1]:
                print('re-initialization was called at iter:' + str(i))
                self.res -= alpha * self._p
                self.initialise_cgls()

            self._r -= alpha * q

            s = Atb(self._r, self.geo, self.angles)
            s_norm = np.linalg.norm(s.ravel(), 2)

            gamma1 = s_norm * s_norm
            beta = gamma1 / self._gamma
            if self.log_parameters:
                self.parameter_history['alpha'][i] = alpha
                self.parameter_history['beta'][i] = beta
                self.parameter_history['gamma'][i] = self._gamma
                self.parameter_history['q_norm'][i] = q_norm
                self.parameter_history['s_norm'][i] = s_norm

            self._gamma = gamma1
            self._p = s + beta * self._p


cgls = decorator(CGLS)

