from __future__ import division

import time
import copy
import numpy as np
import tigre
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.Atb import Atb
from tigre.utilities.Ax import Ax
import tigre.algorithms as algs

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

    def initialize_algo(self):
        self.__r__ = self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        self.__p__ = Atb(self.__r__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
        p_norm = np.linalg.norm(self.__p__.ravel(), 2)
        self.__gamma__ = p_norm * p_norm

    # Overide
    def run_main_iter(self):

        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        avgtime = []

        self.initialize_algo()
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

            self.l2l[0, i] = np.linalg.norm(self.proj - tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids))
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                self.res -= alpha * self.__p__

                if self.verbose:
                    print("re-initilization of CGLS called at iteration:" + str(i))
                if self.re_init_at_iteration + 1 == i or not self.restart:
                    print("CGLS exited due to divergence.")
                    return self.res
                self.re_init_at_iteration=i
                i=i-1
                self.initialize_algo()
                break

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

    def initialize_algo(self):
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
        self.initialize_algo()
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

            self.l2l[0, i] = np.linalg.norm(self.proj - tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids))
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                self.res -= (phi / rho) * (self.__v__-self.__w__)/((theta / rho))
                if self.verbose:
                    print("re-initilization of LSQR called at iteration:" + str(i))
                if self.re_init_at_iteration + 1 == i or not self.restart:
                    print("LSQR exited due to divergence.")
                    return self.res
                self.re_init_at_iteration=i
                i=i-1
                self.initialize_algo()
                break

lsqr = decorator(LSQR, name="lsqr")

class hybrid_LSQR(IterativeReconAlg): 
    __doc__ = (
        " hybrid_LSQR solves the CBCT problem using the hybrid_LSQR\n"
        "  hybrid_LSQR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
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
        
        normr = np.linalg.norm(self.__u__.ravel(), 2)
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

   
            alpha = np.linalg.norm(v.ravel(), 2)
            v = v/alpha
            self.__V__[i] = v.ravel()

            #% Update U_{ii+1}
            self.__u__ = tigre.Ax(v, self.geo, self.angles, "Siddon", gpuids=self.gpuids) - alpha*self.__u__
            
            for j in range(i-1):
                self.__u__=np.reshape(self.__u__.ravel()-(self.__U__[j]*self.__u__.ravel())*self.__U__[j],self.__u__.shape)
                
            
            self.__beta__  = np.linalg.norm(self.__u__.ravel(), 2)
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

      


            self.l2l[0, i] = np.linalg.norm(self.proj - tigre.Ax(self.res + np.reshape(np.matmul(np.transpose(self.__V__[0:i+1]),y),self.res.shape), self.geo, self.angles, "Siddon", gpuids=self.gpuids))
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                if self.re_init_at_iteration + 1 == i or not self.restart:
                    print("hybrid LSQR exited due to divergence at iteration "+str(i))
                    return  self.res + np.reshape(np.matmul(np.transpose(self.__V__[0:i+1]),y),self.res.shape)
                
            #% Test for convergence. 
            #% msl: I still need to implement this. 
            #% msl: There are suggestions on the original paper. Let's talk about it!
        print(y)
        self.res = self.res + np.reshape(np.matmul(np.transpose(self.__V__),y),self.res.shape)
        return self.res
hybrid_lsqr = decorator(hybrid_LSQR, name="hybrid_lsqr")

class LSMR(IterativeReconAlg): 
    __doc__ = (
        " LSMR solves the CBCT problem using LSMR\n"
        "  LSMR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    def initialize_algo(self):
        #% David Chin-Lung Fong and Michael Saunders //doi.org/10.1137/10079687X
        #% Enumeration as given in the paper for 'Algorithm LSMR'
        #% (1) Initialize 
        self.__u__=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        normr = np.linalg.norm(self.__u__.ravel(), 2)
        self.__beta__ = normr
        self.__u__ = self.__u__/normr

        self.__v__ = Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
        self.__alpha__ =np.linalg.norm(self.__v__.ravel(), 2)
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
        self.initialize_algo()
        for i in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(i)
            if self.Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            avgtic = default_timer()  
                    
            #% (3) Continue the bidiagonalization
            self.__u__ = tigre.Ax(self.__v__, self.geo, self.angles, "Siddon", gpuids=self.gpuids) - self.__alpha__*self.__u__
            self.__beta__ = np.linalg.norm(self.__u__.ravel(),2)
            self.__u__ = self.__u__ / self.__beta__
            
            self.__v__ = tigre.Atb(self.__u__, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids) - self.__beta__*self.__v__
            self.__alpha__ = np.linalg.norm(self.__v__.ravel(),2)
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

            self.l2l[0, i] = np.linalg.norm(self.proj - tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids))
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                self.res -= (self.__zeta__ / (self.__rho__*self.__rhobar__)) * self.__hbar__ 
                if self.verbose:
                    print("re-initilization of LSMR called at iteration:" + str(i))
                if self.re_init_at_iteration + 1 == i or not self.restart:
                    print("LSMR exited due to divergence.")
                    return self.res
                self.re_init_at_iteration=iter
                iter=iter-1
                self.initialize_algo()
                break
lsmr = decorator(LSMR, name="lsmr")

class IRN_TV_CGLS(IterativeReconAlg):
    __doc__ = (
        " IRN_TV_CGLS solves the CBCT problem using CGLS with TV constraints\n"
        "  IRN_TV_CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):
        # Don't precompute V and W.
        kwargs.update(dict(W=None, V=None))
        kwargs.update(dict(blocksize=angles.shape[0]))
        self.re_init_at_iteration = 0
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)

    
    def __build_weights__(self):
        Dxx=np.copy(self.res)
        Dyx=np.copy(self.res)
        Dzx=np.copy(self.res)

        Dxx[0:-2,:,:]=self.res[0:-2,:,:]-self.res[1:-1,:,:]
        Dyx[:,0:-2,:]=self.res[:,0:-2,:]-self.res[:,1:-1,:]
        Dzx[:,:,0:-2]=self.res[:,:,0:-2]-self.res[:,:,1:-1]
 
        return (Dxx**2+Dyx**2+Dzx**2+1e-6)**(-1/4)

    def Lx(self,W,img):
        Dxx=np.copy(img)
        Dyx=np.copy(img)
        Dzx=np.copy(img)

        Dxx[0:-2,:,:]=img[0:-2,:,:]-img[1:-1,:,:]
        Dyx[:,0:-2,:]=img[:,0:-2,:]-img[:,1:-1,:]
        Dzx[:,:,0:-2]=img[:,:,0:-2]-img[:,:,1:-1]

        return np.stack((W*Dxx,W*Dyx,W*Dzx),axis=0)
        
    def Ltx(self,W,img3):
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

    def run_main_iter(self):
        self.l2l = np.zeros((1, self.niter*self.niter_outer), dtype=np.float32)
        avgtime = []

        res0=self.res
   
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
            self.res=res0
    
            prox_aux_1 =Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
            prox_aux_2 = self.Lx(W,self.res)*np.sqrt(self.lmbda)
    
            r_aux_1 = self.proj - prox_aux_1
            r_aux_2 = -prox_aux_2
            #% Malena: changed the format, r_aux_2 is 3
            #% r = cat(3,r_aux_1, r_aux_2); % Malena: size guide, erase later, N x N x (100 + N-1)
            p_aux_1 = tigre.Atb(r_aux_1, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
            p_aux_2 = np.sqrt(self.lmbda)*self.Ltx(W, r_aux_2)
            p = p_aux_1 + p_aux_2

            gamma=np.linalg.norm(p.ravel(),2)**2

            for i in range(self.niter):
                res0=self.res

                q_aux_1 = tigre.Ax(p, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
                q_aux_2 = self.Lx(W,p)*np.sqrt(self.lmbda)

                #% q = cat(3, q_aux_1, q_aux_2{1},q_aux_2{2},q_aux_2{3}); % Probably never need to actually do this
                #% alpha=gamma/norm(q(:),2)^2;
                alpha=gamma/(np.linalg.norm(q_aux_1.ravel(),2)**2 + np.linalg.norm(q_aux_2[0].ravel(),2)**2 + np.linalg.norm(q_aux_2[1].ravel(),2)**2+np.linalg.norm(q_aux_2[2].ravel(),2)**2)
                self.res=self.res+alpha*p
                aux=self.proj-tigre.Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
                #% residual norm or the original least squares (not Tikhonov).
                #% Think if that is what we want of the NE residual
                self.l2l[0, outer*self.niter+i] = np.linalg.norm(aux.ravel(),2)

              
                #% If step is adecuate, then continue withg CGLS
                r_aux_1 = r_aux_1-alpha*q_aux_1
                r_aux_2=r_aux_2-alpha*q_aux_2

                s_aux_1 = tigre.Atb(r_aux_1, self.geo, self.angles, backprojection_type="matched", gpuids=self.gpuids)
                s_aux_2 =  np.sqrt(self.lmbda) * self.Ltx(W, r_aux_2)
                s = s_aux_1 + s_aux_2

                gamma1=np.linalg.norm(s.ravel(),2)**2
                beta=gamma1/gamma
                gamma=gamma1
                p=s+beta*p

            avgtoc = default_timer()
            avgtime.append(abs(avgtic - avgtoc))

            if self.Quameasopts is not None:
                self.error_measurement(res_prev, outer)

irn_tv_cgls = decorator(IRN_TV_CGLS, name="irn_tv_cgls")


class AB_GMRES(IterativeReconAlg): 
    __doc__ = (
        " AB_GMRES solves the CBCT problem using preconditioned GMRES\n"
        "  AB_GMRES(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
            # Don't precompute V and W.
            kwargs.update(dict(W=None, V=None))
            kwargs.update(dict(blocksize=angles.shape[0]))
            self.re_init_at_iteration = 0
            IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
            backproject=kwargs.get("backprojector","matched")
            if backproject is "matched":
                self.backproject=Atb
            elif backproject is "FDK":
                self.backproject=algs.fdk

    def __compute_res__(self,x,w,y):
        y=y.astype(np.float32)
        for i in range(w.shape[0]):
            x=x+self.backproject(np.reshape(w[i],self.proj.shape),self.geo,self.angles,gpuids=self.gpuids)*y[i]
        return x

    def run_main_iter(self):

        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        w = np.zeros((self.niter+1,np.prod(self.geo.nDetector)*len(self.angles)),dtype=np.float32)
        r=self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids)
        w[0] = r.ravel()/np.linalg.norm(r.ravel(), 2)
        h=np.zeros((self.niter,self.niter+1),dtype=np.float32)
        for k in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(k)
                
            qk=Ax(self.backproject(np.reshape(w[k],self.proj.shape),self.geo,self.angles,gpuids=self.gpuids),self.geo, self.angles, "Siddon", gpuids=self.gpuids)
            e1=np.zeros(k+2)
            e1[0]=1
            for i in range(k+1):
                h[k,i]=sum(qk.ravel()*w[i])
                qk=qk.ravel()-h[k,i]*w[i]
            
            h[k,k+1]=np.linalg.norm(qk.ravel(),2)
            w[k+1]=qk.ravel()/h[k,k+1]
            y=np.linalg.lstsq(np.transpose(h[0:k+1,0:k+2]),e1*np.linalg.norm(r.ravel(),2),rcond=None)
            y=y[0]
            self.l2l[0, i] = np.linalg.norm((self.proj - tigre.Ax(self.__compute_res__(self.res,w[0:k+1],y),self.geo,self.angles, "Siddon",gpuids=self.gpuids)).ravel(),2)
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                if self.re_init_at_iteration + 1 == i or not self.restart:
                    print("AB-GMRES exited due to divergence at iteration "+str(i))
                    return  self.__compute_res__(self.res,w[0:k+1],y)
             
        self.res=self.__compute_res__(self.res,w[0:-1],y)
        return self.res
ab_gmres = decorator(AB_GMRES, name="ab_gmres")



class BA_GMRES(IterativeReconAlg): 
    __doc__ = (
        " BA_GMRES solves the CBCT problem using preconditioned GMRES\n"
        "  AB_GMRES(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem\n"
        "  using the projection data PROJ taken over ALPHA angles, corresponding\n"
        "  to the geometry descrived in GEO, using NITER iterations."
    ) + IterativeReconAlg.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs): 
            # Don't precompute V and W.
            kwargs.update(dict(W=None, V=None))
            kwargs.update(dict(blocksize=angles.shape[0]))
            self.re_init_at_iteration = 0
            IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
            backproject=kwargs.get("backprojector","matched")
            if backproject is "matched":
                self.backproject=Atb
            elif backproject is "FDK":
                self.backproject=algs.fdk

    def __compute_res__(self,x,w,y):
        y=y.astype(np.float32)
        for i in range(w.shape[0]):
            x=x+np.reshape(w[i],self.res.shape)*y[i]
        return x

    def run_main_iter(self):

        self.l2l = np.zeros((1, self.niter), dtype=np.float32)
        w = np.zeros((self.niter+1,(np.prod(self.geo.nVoxel))),dtype=np.float32)
        r=self.backproject(self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids), self.geo, self.angles, gpuids=self.gpuids)
        w[0] = r.ravel()/np.linalg.norm(r.ravel(), 2)
        h=np.zeros((self.niter,self.niter+1),dtype=np.float32)

        for k in range(self.niter):
            if self.verbose:
                self._estimate_time_until_completion(k)
                
            qk=self.backproject(Ax(np.reshape(w[k],self.res.shape),self.geo,self.angles, "Siddon",gpuids=self.gpuids),self.geo, self.angles, gpuids=self.gpuids)
            e1=np.zeros(k+2)
            e1[0]=1
            for i in range(k+1):
                h[k,i]=sum(qk.ravel()*w[i])
                qk=qk.ravel()-h[k,i]*w[i]
            
            h[k,k+1]=np.linalg.norm(qk.ravel(),2)
            w[k+1]=qk.ravel()/h[k,k+1]
            y=np.linalg.lstsq(np.transpose(h[0:k+1,0:k+2]),e1*np.linalg.norm(r.ravel(),2),rcond=None)
            y=y[0]

            self.l2l[0, i] = np.linalg.norm((self.proj - tigre.Ax(self.__compute_res__(self.res,w[0:k+1],y), self.geo, self.angles, "Siddon", gpuids=self.gpuids)).ravel(),2)
            if i > 0 and self.l2l[0, i] > self.l2l[0, i - 1]:
                if self.re_init_at_iteration + 1 == i or not self.restart:
                    print("BA-GMRES exited due to divergence at iteration "+str(i))
                    return  self.__compute_res__(self.res,w[0:k+1],y)
             
        self.res=self.__compute_res__(self.res,w[0:-1],y)
        return self.res
ba_gmres = decorator(BA_GMRES, name="ba_gmres")