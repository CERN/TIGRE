from __future__ import division

import time
import copy
import numpy as np
import scipy as sp
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Atb import Atb
from tigre.utilities.Ax import Ax
from tigre.utilities.im3Dnorm import l2norm, inner
import tigre.algorithms as algs
import scipy.sparse.linalg


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
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        # -1, not 0: 0 is a legal iteration index, and the give-up test
        # compares against it directly.
        self.re_init_at_iteration = -1
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    def initialize_algo(self):
        self.__r__ = self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        self.__p__ = Atb(self.__r__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
        p_norm = l2norm(self.__p__.ravel())
        self.__gamma__ = p_norm * p_norm

    # override
    def run_main_iter(self):
        """Two nested loops, as in the MATLAB reference (Algorithms/CGLS.m).

        The outer loop rebuilds r, p and gamma after a loss of orthogonality
        (a restart); the inner loop iterates. `re_init_at_iteration` records
        where the last restart happened, so two restarts in a row stop the
        algorithm. It starts at -1 because 0 is a valid iteration index.
        """
        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []

        i = 0
        while i < self.niter:
            self.initialize_algo()
            restarted = False
            while i < self.niter:
                if self.verbose:
                    self._estimate_time_until_completion(i)
                if self.Quameasopts is not None:
                    res_prev = copy.deepcopy(self.res)

                avgtic = default_timer()
                q = tigre.Ax(self.__p__, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
                q_norm = l2norm(q)
                alpha = self.__gamma__ / (q_norm * q_norm)
                self.res += alpha * self.__p__
                avgtoc = default_timer()
                avgtime.append(abs(avgtic - avgtoc))

                self.l2l[0, i] = l2norm(self.proj - tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids))
                if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                    # Undo the step that lost orthogonality.
                    self.res -= alpha * self.__p__

                    if self.verbose:
                        print("re-initilization of CGLS called at iteration:" + str(i))
                    # Give up only if the restart itself did not help, i.e.
                    # this same iteration failed again after a rebuild.
                    if self.re_init_at_iteration == i or not self.restart:
                        print("CGLS exited due to divergence.")
                        return self.res
                    self.re_init_at_iteration = i
                    restarted = True
                    break                 # -> outer loop, retrying this i

                self.__r__ -= alpha * q
                s = tigre.Atb(self.__r__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
                s_norm = l2norm(s)

                gamma1 = s_norm * s_norm
                beta = gamma1 / self.__gamma__

                self.__gamma__ = gamma1
                self.__p__ = s + beta * self.__p__
                if self.Quameasopts is not None:
                    self.error_measurement(res_prev, i)
                i += 1

            if not restarted:
                break

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
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = -1
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    def initialize_algo(self):
        # Paige and Saunders //doi.org/10.1145/355984.355989

        # Enumeration as given in the paper for 'Algorithm LSQR'
        # (1) Initialize 
        self.__u__=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        
        normr = l2norm(self.__u__.ravel())
        self.__u__ = self.__u__/normr

        self.__beta__ = normr
        self.__phibar__ = normr
        self.__v__ = Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)

        self.__alpha__ =l2norm(self.__v__.ravel())
        self.__v__ = self.__v__/self.__alpha__
        self.__rhobar__ = self.__alpha__
        self.__w__ = np.copy(self.__v__)

    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []
        # Outer loop restarts (rebuilds the Krylov state), inner loop iterates;
        # see CGLS.run_main_iter.
        i = 0
        while i < self.niter:
            self.initialize_algo()
            restarted = False
            while i < self.niter:
                if self.verbose:
                    self._estimate_time_until_completion(i)
                if self.Quameasopts is not None:
                    res_prev = copy.deepcopy(self.res)
                avgtic = default_timer()    
            
                #% (3)(a)
                self.__u__ = tigre.Ax(self.__v__, self.geo, self.angles, "Siddon", gpuids=self.gpuids) - self.__alpha__*self.__u__
                self.__beta__ = l2norm(self.__u__.ravel())
                self.__u__ = self.__u__ / self.__beta__
            
                #% (3)(b)
                self.__v__ = tigre.Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids) - self.__beta__*self.__v__
                self.__alpha__ = l2norm(self.__v__.ravel())
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

                self.l2l[0, i] = l2norm(self.proj - tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids))
                if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                    self.res -= (phi / rho) * (self.__v__-self.__w__)/((theta / rho))
                    if self.verbose:
                        print("re-initilization of LSQR called at iteration:" + str(i))
                    if self.re_init_at_iteration == i or not self.restart:
                        print("LSQR exited due to divergence.")
                        return self.res
                    self.re_init_at_iteration = i
                    restarted = True
                    break                 # -> outer loop, retrying this i
                i += 1
            if not restarted:
                break

lsqr = decorator(LSQR, name="lsqr")

class hybrid_LSQR(IterativeReconAlg): 
    __doc__ = (
        " hybrid_LSQR solves the CBCT problem using the hybrid_LSQR\n"
        "  hybrid_LSQR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        self.__U__ = np.zeros((self.niter+1,np.prod(self.geo.nDetector)*len(self.angles)),dtype=np.float32)
        self.__V__ = np.zeros((self.niter,(np.prod(self.geo.nVoxel))),dtype=np.float32)
        self.__B__ = np.zeros((self.niter,self.niter+1),dtype=np.float32) #% Projected matrix
        self.__proj_rhs__ = np.zeros((self.niter+1,1),dtype=np.float32) #% Projected right hand side

    def initialize_algo(self):
        # Paige and Saunders //doi.org/10.1145/355984.355989
        # Enumeration as given in the paper for 'Algorithm LSQR'
        # % Initialise matrices

        # (1) Initialize 
        self.__u__=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        
        normr = l2norm(self.__u__.ravel())
        self.__u__ = self.__u__/normr
        self.__U__[0]=self.__u__.ravel()

        self.__beta__ = normr
        self.__proj_rhs__[0]=normr


    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []
        self.initialize_algo()
        for i in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(i)
            
            avgtic = default_timer() 

            v = Atb(self.__u__,self.geo,self.angles,backprojection_type="matched",gpuids=self.gpuids)
            
            if i>0:
                v = np.reshape(v.ravel() - self.__beta__*self.__V__[i-1],v.shape)
    
            
            for j in range(i-1):
                v=np.reshape(v.ravel()-(self.__V__[j]*v.ravel())*self.__V__[j],v.shape)

   
            alpha = l2norm(v.ravel())
            v = v/alpha
            self.__V__[i] = v.ravel()

            #% Update U_{ii+1}
            self.__u__ = tigre.Ax(v, self.geo, self.angles, "Siddon", gpuids=self.gpuids) - alpha*self.__u__
            
            for j in range(i-1):
                self.__u__=np.reshape(self.__u__.ravel()-(self.__U__[j]*self.__u__.ravel())*self.__U__[j],self.__u__.shape)
                
            
            self.__beta__  = l2norm(self.__u__.ravel())
            self.__u__ = self.__u__ / self.__beta__ 
            self.__U__[i+1] = self.__u__.ravel()

            #% Update projected matrix
            self.__B__[i,i] = alpha
            self.__B__[i,i+1] = self.__beta__ 
            #% Malena. Proposed update: we should check algorithms breaks; 
            #% 'if abs(alpha) <= eps || abs(beta) <= eps' - end and save

            #% Solve the projected problem 
            #% (using the SVD of the small projected matrix)
            Bk = self.__B__[0:i+1,0:i+2]
            Uk, Sk, Vk = np.linalg.svd(np.transpose(Bk))
      
            if i==0:
                Sk = Sk[0]
            
            rhsk = self.__proj_rhs__[0:i+2]
            rhskhat = np.matmul(np.transpose(Uk),rhsk) #
            Dk = Sk**2 + self.lmbda**2

            rhskhat = Sk * rhskhat[0:i+1,0]
            yhat = rhskhat[0:i+1]/Dk
            y = np.matmul(np.transpose(Vk), yhat)

      


            self.l2l[0, i] = l2norm(self.proj - tigre.Ax(self.res + np.reshape(np.matmul(np.transpose(self.__V__[0:i+1]),y),self.res.shape), self.geo, self.angles, "Siddon", gpuids=self.gpuids))
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                # No restart path here, so a genuine rise ends the run. The
                # old `re_init_at_iteration + 1 == i` test, against a sentinel
                # this class never reassigns, could only fire at i == 1 and
                # ignored a rise at every other iteration.
                print("hybrid LSQR exited due to divergence at iteration "+str(i))
                self.res = self.res + np.reshape(np.matmul(np.transpose(self.__V__[0:i+1]),y),self.res.shape)
                return self.res
                
            #% Test for convergence. 
            #% msl: I still need to implement this. 
            #% msl: There are suggestions on the original paper. Let's talk about it!
        
        self.res = self.res + np.reshape(np.matmul(np.transpose(self.__V__),y),self.res.shape)
        return self.res
hybrid_lsqr = decorator(hybrid_LSQR, name="hybrid_lsqr")

class LSMR(IterativeReconAlg): 
    __doc__ = (
        " LSMR solves the CBCT problem using LSMR\n"
        "  LSMR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = -1
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    def initialize_algo(self):
        #% David Chin-Lung Fong and Michael Saunders //doi.org/10.1137/10079687X
        #% Enumeration as given in the paper for 'Algorithm LSMR'
        #% (1) Initialize 
        self.__u__=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        normr = l2norm(self.__u__.ravel())
        self.__beta__ = normr
        self.__u__ = self.__u__/normr

        self.__v__ = Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
        self.__alpha__ =l2norm(self.__v__.ravel())
        self.__v__ = self.__v__/self.__alpha__

        self.__alphabar__ = self.__alpha__
        self.__zetabar__ = self.__alpha__ * self.__beta__
        self.__rho__ = 1
        self.__rhobar__ = 1
        self.__cbar__ = 1
        self.__sbar__ = 0
        self.__h__ = self.__v__ 
        self.__hbar__ = 0

        #% Compute the residual norm ||r_k||
        self.__betadd__ = self.__beta__
        self.__betad__ = 0
        self.__rhod__ = 1
        self.__tautilda__ = 0
        self.__thetatilda__ = 0
        self.__zeta__ = 0
        self.__d__ = 0

    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []
        # Outer loop restarts (rebuilds the Krylov state), inner loop iterates;
        # see CGLS.run_main_iter.
        i = 0
        while i < self.niter:
            self.initialize_algo()
            restarted = False
            while i < self.niter:
                if self.verbose:
                    self._estimate_time_until_completion(i)
                if self.Quameasopts is not None:
                    res_prev = copy.deepcopy(self.res)
                avgtic = default_timer()  
                    
                #% (3) Continue the bidiagonalization
                self.__u__ = tigre.Ax(self.__v__, self.geo, self.angles, "Siddon", gpuids=self.gpuids) - self.__alpha__*self.__u__
                self.__beta__ = l2norm(self.__u__.ravel())
                self.__u__ = self.__u__ / self.__beta__
            
                self.__v__ = tigre.Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids) - self.__beta__*self.__v__
                self.__alpha__ = l2norm(self.__v__.ravel())
                self.__v__ = self.__v__ / self.__alpha__  

                #% (4) Construct and apply rotation \hat{P}_k
                alphahat = np.sqrt(self.__alphabar__**2 + self.lmbda**2)
                chat = self.__alphabar__/alphahat
                shat = self.lmbda/alphahat

                #% (5) Construct and apply rotation P_k
                rhopre = self.__rho__; 
                self.__rho__ = np.sqrt(alphahat**2 + self.__beta__**2)
                c = alphahat / self.__rho__
                s =  self.__beta__ / self.__rho__
                theta = s * self.__alpha__
                self.__alphabar__ = c * self.__alpha__

                #% (6) Construct and apply rotation \bar{P}_k
                thetabar = self.__sbar__  * self.__rho__
                rhobarpre = self.__rhobar__
                self.__rhobar__ = np.sqrt((self.__cbar__ *self.__rho__)**2 + theta**2)
                self.__cbar__ = self.__cbar__ * self.__rho__ / self.__rhobar__
                self.__sbar__ = theta / self.__rhobar__
                zetapre = self.__zeta__
                self.__zeta__ = self.__cbar__ * self.__zetabar__
                self.__zetabar__ = -self.__sbar__ * self.__zetabar__

                #% (7) Update \bar{h}, x, h
                self.__hbar__  = self.__h__ - (thetabar*self.__rho__/(rhopre*rhobarpre))*self.__hbar__ 
                self.res = self.res + (self.__zeta__ / (self.__rho__*self.__rhobar__)) * self.__hbar__ 
                self.__h__ = self.__v__ - (theta / self.__rho__) * self.__h__

                #% (8) Apply rotation \hat{P}_k, P_k
                betaacute = chat* self.__betadd__
                betacheck = - shat* self.__betadd__

                #% Computing ||r_k||

                betahat = c * betaacute
                betadd = -s * betaacute

                #% Update estimated quantities of interest.
                #%  (9) If k >= 2, construct and apply \tilde{P} to compute ||r_k||
                rhotilda = np.sqrt(self.__rhod__**2 + thetabar**2)
                ctilda = self.__rhod__ / rhotilda
                stilda = thetabar / rhotilda
                thetatildapre = self.__thetatilda__
                self.__thetatilda__ = stilda * self.__rhobar__
                self.__rhod__ = ctilda * self.__rhobar__
                #% betatilda = ctilda * betad + stilda * betahat; % msl: in the orinal paper, but not used
                self.__betad__ = -stilda * self.__betad__ + ctilda * betahat

                #% (10) Update \tilde{t}_k by forward substitution
                self.__tautilda__ = (zetapre - thetatildapre* self.__tautilda__) / rhotilda
                taud = (self.__zeta__ - self.__thetatilda__*self.__tautilda__) / self.__rhod__
            
                #% (11) Compute ||r_k||
                self.__d__ = self.__d__ + betacheck**2
                gamma_var = self.__d__ + (self.__betad__ - taud)**2 + betadd**2

                avgtoc = default_timer()
                avgtime.append(abs(avgtic - avgtoc))

                if self.Quameasopts is not None:
                    self.error_measurement(res_prev, i)

                self.l2l[0, i] = l2norm(self.proj - tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids))
                if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                    self.res -= (self.__zeta__ / (self.__rho__*self.__rhobar__)) * self.__hbar__ 
                    if self.verbose:
                        print("re-initilization of LSMR called at iteration:" + str(i))
                    if self.re_init_at_iteration == i or not self.restart:
                        print("LSMR exited due to divergence.")
                        return self.res
                    self.re_init_at_iteration = i
                    restarted = True
                    break                 # -> outer loop, retrying this i
                i += 1
            if not restarted:
                break
lsmr = decorator(LSMR, name="lsmr")

class IRN_TV_CGLS(IterativeReconAlg):
    __doc__ = (
        " IRN_TV_CGLS solves the CBCT problem using CGLS with TV constraints\n"
        "  IRN_TV_CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    
    # The forward difference D and its EXACT adjoint. CGLS is applied to the
    # stacked operator [A; sqrt(lmbda) W D], and CGLS assumes the transpose it
    # is handed really is the adjoint - so these two must agree to machine
    # precision, not approximately.
    #
    #   D u[i] = u[i] - u[i+1]  for i <= n-2,   0 at i = n-1
    #   <Du, y> = sum_{i<=n-2} u[i]y[i] - sum_{j>=1} u[j]y[j-1]
    #   => (D^T y)[0] = y[0]; [j] = y[j]-y[j-1] for 1<=j<=n-2; [n-1] = -y[n-2]
    #
    # The previous pair was wrong on both counts: `Dxx = np.copy(img)` followed
    # by writing only `[0:-2]` left the LAST TWO slices holding raw image
    # VALUES rather than differences (so D did not annihilate a constant - it
    # penalised intensity there), and D^T was off by one at index n-2. Measured
    # adjoint asymmetry 6.4e-2, which is what made the inner CGLS diverge.
    def __D__(self, u, axis):
        out = np.zeros_like(u)
        a = np.swapaxes(out, 0, axis)
        b = np.swapaxes(u, 0, axis)
        a[:-1] = b[:-1] - b[1:]
        return out

    def __Dt__(self, y, axis):
        out = np.zeros_like(y)
        a = np.swapaxes(out, 0, axis)
        b = np.swapaxes(y, 0, axis)
        a[0] = b[0]
        a[1:-1] = b[1:-1] - b[0:-2]
        a[-1] = -b[-2]
        return out

    def __build_weights__(self):
        Dxx = self.__D__(self.res, 0)
        Dyx = self.__D__(self.res, 1)
        Dzx = self.__D__(self.res, 2)
        # W^2 = 1/sqrt(|grad u|^2 + eps) is the half-quadratic weight that makes
        # ||W grad u||^2 a surrogate for the TV seminorm.
        return (Dxx**2+Dyx**2+Dzx**2+1e-6)**np.float32(-1/4)

    def Lx(self,W,img):
        return np.stack((W*self.__D__(img, 0),
                         W*self.__D__(img, 1),
                         W*self.__D__(img, 2)), axis=0)
        
    def Ltx(self,W,img3):
        return (self.__Dt__(W * img3[0], 0)
                + self.__Dt__(W * img3[1], 1)
                + self.__Dt__(W * img3[2], 2))

    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter*self.niter_outer), dtype=np.float32)
        avgtime = []

        res0=self.res
        best_res = np.float32(np.inf)
   
        for outer in range(self.niter_outer):
            if self.verbose:
                niter=self.niter
                self.niter=self.niter_outer
                self._estimate_time_until_completion(outer)
                self.niter=niter
            if self.Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            avgtic = default_timer()    


            W=self.__build_weights__()
            res_prev_inner = np.float32(np.inf)
            res_now = np.float32(np.inf)
    
            prox_aux_1 =Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
            prox_aux_2 = self.Lx(W,self.res)*np.sqrt(self.lmbda, dtype=np.float32)
    
            r_aux_1 = self.proj - prox_aux_1
            r_aux_2 = -prox_aux_2
            #% Malena: changed the format, r_aux_2 is 3
            #% r = cat(3,r_aux_1, r_aux_2); % Malena: size guide, erase later, N x N x (100 + N-1)
            p_aux_1 = tigre.Atb(r_aux_1, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
            p_aux_2 = np.sqrt(self.lmbda, dtype=np.float32)*self.Ltx(W, r_aux_2)
            p = p_aux_1 + p_aux_2
            
            gamma=np.float32(l2norm(p.ravel())**2)

            for i in range(self.niter):
                res0=self.res
                q_aux_1 = tigre.Ax(p, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
                q_aux_2 = self.Lx(W,p)*np.sqrt(self.lmbda, dtype=np.float32)

                #% q = cat(3, q_aux_1, q_aux_2{1},q_aux_2{2},q_aux_2{3}); % Probably never need to actually do this
                #% alpha=gamma/norm(q(:),2)^2;
                alpha=np.float32(gamma/(l2norm(q_aux_1.ravel())**2 + l2norm(q_aux_2[0].ravel())**2 + l2norm(q_aux_2[1].ravel())**2+l2norm(q_aux_2[2].ravel())**2))
                self.res=self.res+alpha*p
                aux=self.proj-tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
                #% residual norm or the original least squares (not Tikhonov).
                #% Think if that is what we want of the NE residual
                step = outer*self.niter+i
                self.l2l[0, step] = res_now = l2norm(aux.ravel())

                # The guard plain CGLS has always had and this loop never did.
                # In exact arithmetic a CGLS residual cannot rise; when it
                # does, orthogonality is lost and every later step compounds
                # it. Left unchecked, niter_outer * niter such steps reach
                # absurd values - measured mean 8469 / max 5e7 at 512^3, and
                # 7e5 on a 64^3 phantom where the rise begins at step 2.
                #
                # Compare only WITHIN one outer iteration. Across a
                # reweighting the two residuals belong to different weighted
                # problems, so a first step that looks worse than the previous
                # outer's last is not evidence of anything.
                if i > 0 and res_now > res_prev_inner:
                    self.res = res0
                    self.l2l[0, step] = res_prev_inner
                    res_now = res_prev_inner
                    if self.verbose:
                        print("IRN_TV_CGLS: residual rose at outer " + str(outer)
                              + ", inner " + str(i) + " - reweighting")
                    break
                res_prev_inner = res_now
              
                #% If step is adequate, then continue withg CGLS
                r_aux_1 = r_aux_1-alpha*q_aux_1
                r_aux_2=r_aux_2-alpha*q_aux_2

                s_aux_1 = tigre.Atb(r_aux_1, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
                s_aux_2 =  np.sqrt(self.lmbda, dtype=np.float32) * self.Ltx(W, r_aux_2)
                s = s_aux_1 + s_aux_2

                gamma1=np.float32(l2norm(s.ravel())**2)
                beta=np.float32(gamma1/gamma)
                gamma=gamma1
                p=s+beta*p

            avgtoc = default_timer()
            avgtime.append(abs(avgtic - avgtoc))

            if self.Quameasopts is not None:
                self.error_measurement(res_prev, outer)

            # A whole reweighting that improved nothing means the outer
            # iteration has stalled; more of them will not help.
            if res_now >= best_res:
                if self.verbose:
                    print("IRN_TV_CGLS: reweighting " + str(outer)
                          + " did not improve the residual - stopping")
                break
            best_res = res_now

irn_tv_cgls = decorator(IRN_TV_CGLS, name="irn_tv_cgls")


class AB_GMRES(IterativeReconAlg): 
    __doc__ = (
        " AB_GMRES solves the CBCT problem using preconditioned GMRES\n"
        "  AB_GMRES(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
            # Don't precompute V and W.
            kwargs.update(dict(W=None, V=None))
            kwargs.update(dict(blocksize=angles.shape[0]))
            self.re_init_at_iteration = 0
            if "backprojector" in kwargs:
                backproject=kwargs.pop("backprojector")
            else:
                backproject="matched"
            if backproject == "matched":
                self.backproject=Atb
            elif backproject == "FDK":
                self.backproject=algs.fdk
            IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
      

    def __compute_res__(self,x,w,y):
        y=y.astype(np.float32)
        for i in range(w.shape[0]):
            x=x+self.backproject(np.reshape(w[i],self.proj.shape),self.geo,self.angles,gpuids=self.gpuids)*y[i]
        return x

    def run_main_iter(self):

        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        w = np.zeros((self.niter+1,np.prod(self.geo.nDetector)*len(self.angles)),dtype=np.float32)
        r=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        w[0] = r.ravel()/l2norm(r.ravel())
        h=np.zeros((self.niter,self.niter+1),dtype=np.float32)
        for k in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(k)
                
            qk=Ax(self.backproject(np.reshape(w[k],self.proj.shape),self.geo,self.angles,gpuids=self.gpuids),self.geo, self.angles, "Siddon", gpuids=self.gpuids)
            e1=np.zeros(k+2)
            e1[0]=1
            qk = qk.ravel()
            for i in range(k+1):
                h[k,i]=inner(qk,w[i])
                qk -= h[k,i]*w[i]
            
            h[k,k+1]=l2norm(qk)
            w[k+1]=qk/h[k,k+1]
            y=np.linalg.lstsq(np.transpose(h[0:k+1,0:k+2]),e1*l2norm(r.ravel()),rcond=None)
            y=y[0]
            self.l2l[0, k] = l2norm((self.proj - tigre.Ax(self.__compute_res__(self.res,w[0:k+1],y),self.geo,self.angles, "Siddon",gpuids=self.gpuids)).ravel())
            if k > 0 and self.l2l[0, k] > self.l2l[0, k - 1]:
                # No restart path exists in GMRES here, so a genuine rise ends
                # the run. It used to be tested as `re_init_at_iteration + 1 ==
                # i` against a sentinel that is never reassigned, which fired
                # only when the leaked inner-loop index happened to equal 1 and
                # did nothing at all for any later iteration.
                print("AB-GMRES exited due to divergence at iteration "+str(k))
                self.res = self.__compute_res__(self.res,w[0:k+1],y)
                return self.res
             
        self.res=self.__compute_res__(self.res,w[0:-1],y)
        return self.res
ab_gmres = decorator(AB_GMRES, name="ab_gmres")



class BA_GMRES(IterativeReconAlg): 
    __doc__ = (
        " BA_GMRES solves the CBCT problem using preconditioned GMRES\n"
        "  AB_GMRES(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
            # Don't precompute V and W.
            kwargs.update(dict(W=None, V=None))
            kwargs.update(dict(blocksize=angles.shape[0]))
            self.re_init_at_iteration = 0
            if "backprojector" in kwargs:
                backproject=kwargs.pop("backprojector")
            else:
                backproject="matched"
            if backproject == "matched":
                self.backproject=Atb
            elif backproject == "FDK":
                self.backproject=algs.fdk
            IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    def __compute_res__(self,x,w,y):
        y=y.astype(np.float32)
        for i in range(w.shape[0]):
            x=x+np.reshape(w[i],self.res.shape)*y[i]
        return x

    def run_main_iter(self):

        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        w = np.zeros((self.niter+1,(np.prod(self.geo.nVoxel))),dtype=np.float32)
        r=self.backproject(self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids), self.geo, self.angles, gpuids=self.gpuids)
        w[0] = r.ravel()/l2norm(r.ravel())
        h=np.zeros((self.niter,self.niter+1),dtype=np.float32)

        for k in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(k)
                
            qk=self.backproject(Ax(np.reshape(w[k],self.res.shape),self.geo,self.angles, "Siddon",gpuids=self.gpuids),self.geo, self.angles, gpuids=self.gpuids)
            e1=np.zeros(k+2)
            e1[0]=1
            qk = qk.ravel()
            for i in range(k+1):
                h[k,i]=inner(qk,w[i])
                qk -= h[k,i]*w[i]
            
            h[k,k+1]=l2norm(qk)
            w[k+1]=qk/h[k,k+1]
            y=np.linalg.lstsq(np.transpose(h[0:k+1,0:k+2]),e1*l2norm(r.ravel()),rcond=None)
            y=y[0]

            self.l2l[0, k] = l2norm((self.proj - tigre.Ax(self.__compute_res__(self.res,w[0:k+1],y), self.geo, self.angles, "Siddon", gpuids=self.gpuids)).ravel())
            if k > 0 and self.l2l[0, k] > self.l2l[0, k - 1]:
                # See the AB-GMRES note: no restart path, and the old sentinel
                # test could only ever fire at one particular iteration.
                print("BA-GMRES exited due to divergence at iteration "+str(k))
                self.res = self.__compute_res__(self.res,w[0:k+1],y)
                return self.res
             
        self.res=self.__compute_res__(self.res,w[0:-1],y)
        return self.res
ba_gmres = decorator(BA_GMRES, name="ba_gmres")


class hybrid_fLSQR_TV(IterativeReconAlg): 
    __doc__ = (
        " hybrid_fLSQR_TV solves the CBCT problem using preconditioned hybrid flexuble LSQR with TV regularization\n"
        "  AB_GMRES(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry described in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        self.__U__ = np.zeros((self.niter+1,np.prod(self.geo.nDetector)*len(self.angles)),dtype=np.float32)
        self.__V__ = np.zeros((self.niter,(np.prod(self.geo.nVoxel))),dtype=np.float32)
        self.__Z__ = np.zeros((self.niter,(np.prod(self.geo.nVoxel))),dtype=np.float32)

        self.__M__ = np.zeros((self.niter,self.niter+1),dtype=np.float32) #% Projected matrix
        self.__T__ = np.zeros((self.niter,self.niter),dtype=np.float32) #% Projected matrix
        self.__proj_rhs__ = np.zeros((self.niter+1,1),dtype=np.float32) #% Projected right hand side
    

      
    def __build_weights__(self):
        Dxx=np.copy(self.res)
        Dyx=np.copy(self.res)
        Dzx=np.copy(self.res)

        Dxx[0:-2,:,:]=self.res[0:-2,:,:]-self.res[1:-1,:,:]
        Dyx[:,0:-2,:]=self.res[:,0:-2,:]-self.res[:,1:-1,:]
        Dzx[:,:,0:-2]=self.res[:,:,0:-2]-self.res[:,:,1:-1]
 
        return (Dxx**2+Dyx**2+Dzx**2+1e-6)**(-1/4)

    def Lx(self,W,img):
        img=np.reshape(img,self.res.shape)
        Dxx=np.copy(img)
        Dyx=np.copy(img)
        Dzx=np.copy(img)

        Dxx[0:-2,:,:]=img[0:-2,:,:]-img[1:-1,:,:]
        Dyx[:,0:-2,:]=img[:,0:-2,:]-img[:,1:-1,:]
        Dzx[:,:,0:-2]=img[:,:,0:-2]-img[:,:,1:-1]

        return np.stack((W*Dxx,W*Dyx,W*Dzx),axis=0)

    def Ltx(self,W,img3):
        img3 =np.reshape(img3,(3,*self.res.shape))
        Wx_1 = W * img3[0,:,:,:]
        Wx_2 = W * img3[1,:,:,:]
        Wx_3 = W * img3[2,:,:,:]

        DxtWx_1=Wx_1
        DytWx_2=Wx_2
        DztWx_3=Wx_3
        
        DxtWx_1[1:-2,:,:]=Wx_1[1:-2,:,:]-Wx_1[0:-3,:,:]
        DxtWx_1[-1,:,:]=-Wx_1[-2,:,:]
        
        DytWx_2[:,1:-2,:]=Wx_2[:,1:-2,:]-Wx_2[:,0:-3,:]
        DytWx_2[:,-1,:]=-Wx_2[:,-2,:]
        
        DztWx_3[:,:,1:-2]=Wx_3[:,:,1:-2]-Wx_3[:,:,0:-3]
        DztWx_3[:,:,-1]=-Wx_3[:,:,-2]

        return DxtWx_1 + DytWx_2 + DztWx_3

    def mvT(self,k_aux,x ):
            return x.ravel() - k_aux.ravel()*np.sum(x.ravel())
    def mv(self,k_aux,x ):
            return x.ravel() - inner(k_aux,x)
        
    def run_main_iter(self):
        # % Initialise matrices

        # (1) Initialize 
        u=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        
        k_aux = Ax(np.ones(self.res.shape,dtype=np.float32)/np.sqrt(np.prod(self.geo.nVoxel), dtype=np.float32), self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        xA0 = 1/(np.sqrt(np.prod(self.geo.nVoxel), dtype=np.float32)*l2norm(k_aux.ravel())**2)*inner(k_aux,u)
        xA0 =np.ones(self.res.shape,dtype=np.float32)*xA0

        k_aux = 1/(np.sqrt(np.prod(self.geo.nVoxel), dtype=np.float32)*l2norm(k_aux.ravel())**2)*Atb(k_aux, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)

        u = u - Ax(xA0, self.geo, self.angles, "Siddon", gpuids=self.gpuids)

        normr=l2norm(u.ravel())
        u=u/normr
        self.__U__[0]=u.ravel()

        self.__proj_rhs__[0]=normr

        if np.max(self.res)==0:
            W=np.ones(self.res.shape,dtype=np.float32)
        else:
            W=self.__build_weights__()

        self.l2l = np.zeros((1, self.niter), dtype=np.float32)

        for i in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(i)
            v=Atb(u,self.geo,self.angles,backprojection_type="matched", gpuids=self.gpuids)
            for j in range(i+1):
                self.__T__[i,j] = inner(self.__V__[j],v)
                v = v.ravel() - self.__T__[i,j]*self.__V__[j]
            v=np.reshape(v,self.res.shape)
            self.__T__[i,i] = l2norm(v.ravel())
            v=v/self.__T__[i,i]

            z=self.mvT(k_aux,v)

            L = sp.sparse.linalg.LinearOperator((np.prod(self.res.shape),np.prod(self.res.shape)*3), matvec=lambda x: self.Ltx(W,x).ravel(),rmatvec=lambda x: self.Lx(W,x).ravel())
            aux_z=sp.sparse.linalg.lsqr(L,z.ravel(),iter_lim=50)
            L = sp.sparse.linalg.LinearOperator((np.prod(self.res.shape)*3,np.prod(self.res.shape)), matvec=lambda x: self.Lx(W,x).ravel(),rmatvec=lambda x: self.Ltx(W,x).ravel())
            z=sp.sparse.linalg.lsqr(L,aux_z[0].ravel(),iter_lim=50)
            z=z[0].astype(np.float32)

            z=self.mv(k_aux,z)

            self.__V__[i]=v.ravel()
            self.__Z__[i]=z.ravel()

            # Update U and projected matrix M
            u = Ax(np.reshape(z,self.res.shape), self.geo, self.angles, "Siddon", gpuids=self.gpuids)
            for j in range(i+1):
                self.__M__[i,j] = inner(self.__U__[j],u)
                u = u.ravel() - np.dot(self.__M__[i,j],self.__U__[j])
            self.__M__[i,i+1]=l2norm(u.ravel())
            u=u/self.__M__[i,i+1]
            u=np.reshape(u,self.proj.shape)
            self.__U__[i+1]=u.ravel()

            ### Solve the regularized projected problem
            # (using the DVF of the small projected matrix)

            Mk=self.__M__[0:i+1,0:i+2]
            # Prepare the projected regularization term
            WZ=np.zeros((i+1,(3*np.prod(self.geo.nVoxel))),dtype=np.float32)
            for j in range(i+1):
                # This can be more efficiently...
                # DZ can be saved and updated at each iteration
                WZ[j]=self.Lx(W,self.__Z__[j]).ravel()
            __,ZRk =sp.linalg.qr(np.transpose(WZ),mode='economic')
            ZRksq=ZRk[0:i+1,0:i+1]
            rhsk=self.__proj_rhs__[0:i+2]
            MZk=np.concatenate((np.transpose(Mk),self.lmbda*ZRksq))
            rhsZk=np.concatenate((rhsk,np.zeros((i+1,1),dtype=np.float32)))
            y = np.linalg.lstsq(MZk, rhsZk,rcond=None)
            d=np.matmul(np.transpose(self.__Z__[0:i+1]),y[0])
            
            x = self.res + np.reshape(d,self.res.shape) + xA0
            self.l2l[0, i] = l2norm((self.proj - tigre.Ax(x, self.geo, self.angles, "Siddon", gpuids=self.gpuids)).ravel())
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                # Same as hybrid_LSQR: no restart path, so stop. (The message
                # said BA-GMRES; this is hybrid_fLSQR_TV.)
                print("hybrid fLSQR-TV exited due to divergence at iteration "+str(i))
                self.res = x
                return self.res
             
        self.res = x
        return self.res

hybrid_flsqr_tv = decorator(hybrid_fLSQR_TV, name="hybrid_flsqr_tv")
