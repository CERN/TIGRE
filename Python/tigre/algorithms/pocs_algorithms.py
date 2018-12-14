from __future__ import division
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.utilities.im3Dnorm import im3DNORM
from tigre.algorithms.single_pass_algorithms import FDK
from tigre.utilities.Ax import Ax
import time
import copy
import numpy as np
class ASD_POCS(IterativeReconAlg):
    __doc__ = (" solves the reconstruction problem\n"
               " using the projection data PROJ taken over ALPHA angles, corresponding\n"
               " to the geometry descrived in GEO, using NITER iterations.\n") + IterativeReconAlg.__doc__
    def __init__(self,proj,geo,angles,niter, **kwargs):
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        if not kwargs.has_key('alpha'):
            self.alpha = 0.002
        if not kwargs.has_key('alpha_red'):
            self.alpha_red = 0.95
        if not kwargs.has_key('rmax'):
            self.rmax = 0.95
        if not kwargs.has_key('maxl2err'):
            self.epsilon = im3DNORM(FDK(proj,geo,angles),2)*0.2
        if not kwargs.has_key("numiter_tv"):
            self.numiter_tv = 20
        if not kwargs.has_key('regularisation'):
            self.regularisation = 'minimizeTV'
        self.beta =self.lmbda
        self.beta_red = self.lmbda_red


    #Overide
    def run_main_iter(self):
        stop_criteria = False
        iter = 0
        while not stop_criteria:
            if self.verbose:
                if iter == 0:
                    print("POCS Algorithm in progress.")
                    toc = time.clock()
                if iter == 1:
                    tic = time.clock()
                    print('Esitmated time until completetion (s): ' + str((self.niter - 1) * (tic - toc)))
            res_prev = copy.deepcopy(self.res)
            iter +=1
            getattr(self, self.dataminimizing)()
            g = Ax(self.res,self.geo,self.angles)
            dd = im3DNORM(g-self.proj,2)
            dp_vec = self.res -res_prev
            dp = im3DNORM(dp_vec,2)
            if iter==1:
                dtvg = self.alpha * dp
            res_prev = copy.deepcopy(self.res)
            self.res = getattr(self,self.regularisation)(self.res,dtvg)
            dg_vec = self.res - res_prev
            dg = im3DNORM(dg_vec,2)
            if dg>self.rmax*dp and dd>self.epsilon:
                dtvg = dtvg*self.alpha_red
            self.beta *= self.beta_red
            c = np.dot(dg_vec.reshape(-1,),dp_vec.reshape(-1,))/max((im3DNORM(dg_vec,2)*im3DNORM(dp_vec,2)),1e-6)
            if (c<-0.99 and dd<=self.epsilon) or self.beta<0.005 or iter>self.niter:
                if self.verbose:
                    print("\n"
                          "     Stop criteria met: \n"
                          "     c = "+ str(c) + "\n"
                          "     beta = " + str(self.beta)+"\n" 
                          "     iter = " + str(iter)) +"\n"
                stop_criteria = True



asd_pocs = decorator(ASD_POCS,name='asd_pocs')

class AwASD_POCS(ASD_POCS):
    __doc__ = ASD_POCS.__doc__

    def __init__(self,proj,geo,angles,niter,**kwargs):
        kwargs.update(dict(regularisation = 'minimizeAwTV'))
        if not kwargs.has_key('delta'):
            self.delta = np.array([-0.005],dtype=np.float32)[0]
        ASD_POCS.__init__(self,proj,geo,angles,niter,**kwargs)

awasd_pocs = decorator(AwASD_POCS,name='awasd_pocs')
