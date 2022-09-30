from __future__ import division

import time
import copy
import numpy as np
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Atb import Atb
from tigre.utilities.Ax import Ax


if hasattr(time, "perf_counter"):
    default_timer = time.perf_counter
else:
    default_timer = time.clock


class CGLS(IterativeReconAlg):  # noqa: D101
    __doc__ = (
        " CGLS solves the CBCT problem using the conjugate gradient least\n"
        " squares\n"
        " \n"
        "  CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

        self.__r__ = self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        self.__p__ = Atb(self.__r__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    def reinitialise_cgls(self):
        self.__r__ = self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        self.__p__ = Atb(self.__r__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    # Overide
    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []
        for i in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(i)
            if self.Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)

            avgtic = default_timer()
            q = tigre.Ax(self.__p__, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
            q_norm = np.linalg.norm(q)
            alpha = self.__gamma__ / (q_norm * q_norm)
            self.res += alpha * self.__p__
            avgtoc = default_timer()
            avgtime.append(abs(avgtic - avgtoc))
            for item in self.__dict__:
                if (
                    isinstance(getattr(self, item), np.ndarray)
                    and np.isnan(getattr(self, item)).any()
                ):
                    raise ValueError("nan found for " + item + " at iteraton " + str(i))

            aux = self.proj - tigre.Ax(
                self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids
            )
            self.l2l[0, i] = np.linalg.norm(aux)
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                if self.verbose:
                    print("re-initilization of CGLS called at iteration:" + str(i))
                if self.re_init_at_iteration + 1 == i:
                    if self.verbose:
                        print("Algorithm exited with two consecutive reinitializations.")
                    return self.res
                self.res -= alpha * self.__p__
                self.reinitialise_cgls()
                self.re_init_at_iteration = i

            self.__r__ -= alpha * q
            s = tigre.Atb(self.__r__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
            s_norm = np.linalg.norm(s)

            gamma1 = s_norm * s_norm
            beta = gamma1 / self.__gamma__

            self.__gamma__ = gamma1
            self.__p__ = s + beta * self.__p__
            if self.Quameasopts is not None:
                self.error_measurement(res_prev, i)

        if self.verbose:
            print(
                "Average time taken for each iteration for CGLS:"
                + str(sum(avgtime) / len(avgtime))
                + "(s)"
            )

cgls = decorator(CGLS, name="cgls")

class LSQR(IterativeReconAlg): 
    __doc__ = (
        " LSQR solves the CBCT problem using the  least squares\n"
        "  LSQR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        # Paige and Saunders //doi.org/10.1145/355984.355989

        # Enumeration as given in the paper for 'Algorithm LSQR'
        # (1) Initialize 
        self.__u__=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        
        normr = np.linalg.norm(self.__u__.ravel(), 2)
        self.__u__ = self.__u__/normr

        self.__beta__ = normr
        self.__phibar__ = normr
        self.__v__ = Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)

        self.__alpha__ =np.linalg.norm(self.__v__.ravel(), 2)
        self.__v__ = self.__v__/self.__alpha__
        self.__rhobar__ = self.__alpha__
        self.__w__ = np.copy(self.__v__)

    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []
        for i in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(i)
            if self.Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            avgtic = default_timer()    
            
            #% (3)(a)
            self.__u__ = tigre.Ax(self.__v__, self.geo, self.angles, "Siddon", gpuids=self.gpuids) - self.__alpha__*self.__u__
            self.__beta__ = np.linalg.norm(self.__u__.ravel(),2)
            self.__u__ = self.__u__ / self.__beta__
            
            #% (3)(b)
            self.__v__ = tigre.Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids) - self.__beta__*self.__v__
            self.__alpha__ = np.linalg.norm(self.__v__.ravel(),2)
            self.__v__ = self.__v__ / self.__alpha__    

            #% (4)(a-g)
            rho = np.sqrt(self.__rhobar__**2 + self.__beta__**2)
            c = self.__rhobar__ / rho
            s =  self.__beta__ / rho
            theta = s * self.__alpha__    
            self.__rhobar__ = - c * self.__alpha__    
            phi = c * self.__phibar__
            self.__phibar__ = s * self.__phibar__
            
            #% (5) Update x, w
            self.res = self.res + (phi / rho) * self.__w__
            self.__w__ = self.__v__ - (theta / rho) * self.__w__

            avgtoc = default_timer()
            avgtime.append(abs(avgtic - avgtoc))
            
            if self.Quameasopts is not None:
                self.error_measurement(res_prev, i)
           
lsqr = decorator(LSQR, name="lsqr")
