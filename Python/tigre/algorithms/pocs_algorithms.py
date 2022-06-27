from __future__ import division

import copy

import numpy as np
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from tigre.algorithms.iterative_recon_alg import decorator
from tigre.algorithms.single_pass_algorithms import FDK
from tigre.utilities.Ax import Ax
from tigre.utilities.im3Dnorm import im3DNORM




class initASD_POCS(IterativeReconAlg):  
    """
    ASD_POCS solves the ASD_POCS total variation constrained image in 3D
    tomography

    Parameters
    ----------
    :param proj: (np.ndarray, dtype=np.float32)
    Input data, shape = (geo.nDector, nangles)

    :param geo: (tigre.geometry)
    Geometry of detector and image (see examples/Demo code)

    :param angles: (np.ndarray , dtype=np.float32)
    angles of projection, shape = (nangles,3)

    :param niter: (int)
    number of iterations for reconstruction algorithm

    :param kwargs: (dict)
    optional parameters

    Keyword Arguments
    -----------------
    :keyword blocksize: (int)
        number of angles to be included in each iteration
        of proj and backproj in OS_SART update

    :keyword lmbda: (np.float64)
        Sets the value of the hyperparameter for SART iteration.

    :keyword lmbda_red: (np.float64)
        Reduction of lmbda every iteration
        lmbda=lmbda_red*lmbda. Default is 0.99

    :keyword init: (str)
        Describes different initialization techniques.
            None      : Initializes the image to zeros (default)
            "FDK"      : intializes image to FDK reconstrucition

    :keyword verbose:  (Boolean)
        Feedback print statements for algorithm progress
        default=True

    :keyword OrderStrategy : (str)
        Chooses the subset ordering strategy. Options are:
            "ordered"        : uses them in the input order, but
                            divided
            "random"         : orders them randomply
            
    :keyword tviter: (int)
        Number of iterations of minmizeTV to be performed for each iter-
        ation. Default: 20

    :keyword alpha: (float)
        Defines the TV hyperparameter.
        Default: 0.002

    :keyword alpha_red: (float)
        Defines the reduction rate of the hyperparameter
        Default: 0.95

    :keyword rmax: (float)
        The maximum allowed image/TV update ration. If the TV
        update changes the image more than this, the parameter
        will be reduced. Default: 0.95

    :keyword maxl2err: (float)
        Maximum L2 error to accept an image as valid. This
        parameter is crucial for the algorithm, determines at
        what point an image should not be updated further.
        Default is 20% of the FDK L2 norm.
    --------------------------------------------------------------------
    This file is part of the TIGRE Toolbox

    Copyright (c) 2015, University of Bath and
                        CERN-European Organization for Nuclear Research
                        All rights reserved.

    License:            Open Source under BSD.
                        See the full license at
                        https://github.com/CERN/TIGRE/license.txt

    Contact:            tigre.toolbox@gmail.com
    Codes:              https://github.com/CERN/TIGRE/
    --------------------------------------------------------------------
    Coded by:          MATLAB (original code): Ander Biguri
                        PYTHON : Reuben Lindroos
    """

    def __init__(self, proj, geo, angles, niter, **kwargs):

        
        # kwargs.update(dict(regularisation="minimizeTV"))
        # if "blocksize" not in kwargs:
        #     kwargs.update(dict(blocksize=1))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)       
        self.alpha = 0.002 if "alpha" not in kwargs else kwargs["alpha"]
        self.alpha_red = 0.95 if "alpha_red" not in kwargs else kwargs["alpha_red"]
        self.rmax = 0.95 if "rmax" not in kwargs else kwargs["rmax"]
        if "maxl2err" not in kwargs:
            self.epsilon = (
                im3DNORM(Ax(FDK(proj, geo, angles, gpuids=self.gpuids), geo, angles) - proj, 2)
                * 0.2
            )
        else:
            self.epsilon = kwargs["maxl2err"]
        self.numiter_tv = 20 if "tviter" not in kwargs else kwargs["tviter"]
        self.beta = self.lmbda
        self.beta_red = self.lmbda_red
        

    # Override
    def run_main_iter(self):
        stop_criteria = False
        n_iter = 0
        while not stop_criteria:
            if self.verbose:
                self._estimate_time_until_completion(n_iter)

            res_prev = copy.deepcopy(self.res)
            getattr(self, self.dataminimizing)()
            if self.Quameasopts is not None:
                self.error_measurement(res_prev, n_iter)
            n_iter += 1
            g = Ax(self.res, self.geo, self.angles, gpuids=self.gpuids)
            dd = im3DNORM(g - self.proj, 2)
            dp_vec = self.res - res_prev
            dp = im3DNORM(dp_vec, 2)

            if n_iter == 1:
                dtvg = self.alpha * dp

            res_prev = copy.deepcopy(self.res)
            self.res = getattr(self, self.regularisation)(self.res, dtvg)
            dg_vec = self.res - res_prev
            dg = im3DNORM(dg_vec, 2)

            if dg > self.rmax * dp and dd > self.epsilon:
                dtvg = dtvg * self.alpha_red

            self.beta *= self.beta_red
            c = np.dot(dg_vec.reshape(-1,), dp_vec.reshape(-1,)) / max(
                dg * dp, 1e-6
            )  # reshape ensures no copy is made.
            if (c < -0.99 and dd <= self.epsilon) or self.beta < 0.005 or n_iter >= self.niter:
                if self.verbose:
                    print(
                        "\n"
                        "     Stop criteria met: \n"
                        "     c = " + str(c) + "\n"
                        "     beta = " + str(self.beta) + "\n"
                        "     iter = " + str(n_iter) + "\n"
                    )
                stop_criteria = True


class ASD_POCS(initASD_POCS):
    __doc__ = initASD_POCS.__doc__
    
    def __init__(self, proj, geo, angles, niter, **kwargs):
        
        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to 1, please use an OS version of the algorithm for blocksize > 1')
        kwargs.update(dict(blocksize=1))
        kwargs.update(dict(regularisation="minimizeTV"))
        
        initASD_POCS.__init__(self, proj, geo, angles, niter, **kwargs)
       

asd_pocs = decorator(ASD_POCS, name="asd_pocs")


class AwASD_POCS(initASD_POCS):  
    __doc__ = ("    AwASD_POCS is the Adaptive Weighted TV (edge preserving) version of ASD_POCS\n\n"
               "        :extra kwargs delta: (float)\n"
               "            Control amount of smoothing at edges of the image\n"
               "            Default delta = -0.005\n") + initASD_POCS.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to 1, please use an OS version of the algorithm for blocksize > 1')
        kwargs.update(dict(blocksize=1))
        kwargs.update(dict(regularisation="minimizeAwTV"))
        self.delta = np.float32(-0.005) if "delta" not in kwargs else kwargs["delta"] 
        
        initASD_POCS.__init__(self, proj, geo, angles, niter, **kwargs)


awasd_pocs = decorator(AwASD_POCS, name="awasd_pocs")


class OS_ASD_POCS(initASD_POCS):
    __doc__ = ("    Oriented Subsets version of ASD_POCS \n\n"
               "    Default blocksize = 20\n") + initASD_POCS.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        self.blocksize = 20 if "blocksize" not in kwargs else kwargs["blocksize"]
        kwargs.update(dict(regularisation="minimizeTV"))
        
        initASD_POCS.__init__(self, proj, geo, angles, niter, **kwargs)


os_asd_pocs = decorator(OS_ASD_POCS, name="os_asd_pocs")


class OS_AwASD_POCS(initASD_POCS): 
    __doc__ = (
        "    Oriented Subsets and Adaptive Weighted TV version of AwASD_POCS\n\n"
        "    Default blocksize = 20\n\n"
        "        :extra kwargs delta: (float)\n"
        "            Control amount of smoothing at edges of the image\n"
        "            Default delta = -0.005\n"
        ) + initASD_POCS.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        self.blocksize = 20 if "blocksize" not in kwargs else kwargs["blocksize"]
        kwargs.update(dict(regularisation="minimizeAwTV"))
        self.delta = np.float32(-0.005) if "delta" not in kwargs else kwargs["delta"]
        
        initASD_POCS.__init__(self, proj, geo, angles, niter, **kwargs) 


os_awasd_pocs = decorator(OS_AwASD_POCS, name="os_awasd_pocs")




class initPCSD(IterativeReconAlg):
    """
    PCSD solves the reconstruction problem using projection-controlled 
    steepest descent method

    Parameters
    ----------
    :param proj: (np.ndarray, dtype=np.float32)
    Input data, shape = (geo.nDector, nangles)

    :param geo: (tigre.geometry)
    Geometry of detector and image (see examples/Demo code)

    :param angles: (np.ndarray , dtype=np.float32)
    angles of projection, shape = (nangles,3)

    :param niter: (int)
    number of iterations for reconstruction algorithm

    :param kwargs: (dict)
    optional parameters

    Keyword Arguments
    -----------------
    :keyword blocksize: (int)
        number of angles to be included in each iteration
        of proj and backproj in OS_SART update
        
    :keyword lambda: (np.float64)
        Sets the value of the hyperparameter for SART iteration.

    :keyword lambda_red: (np.float64)
        Reduction of lambda every iteration
        lambda=lambda_red*lambda. Default is 0.99

    :keyword init: (str)
        Describes different initialization techniques.
               None      : Initializes the image to zeros (default)
              "FDK"      : intializes image to FDK reconstrucition

    :keyword verbose:  (Boolean)
        Feedback print statements for algorithm progress
        default=True

    :keyword OrderStrategy : (str)
        Chooses the subset ordering strategy. Options are:
                 "ordered"        : uses them in the input order, but
                                    divided
                 "random"         : orders them randomply
                 
    :keyword tviter: (int)
        Number of iterations of minmizeTV to be performed for each iter-
        ation. Default: 20

    :keyword maxl2err: (float)
         Maximum L2 error to accept an image as valid. This
         parameter is crucial for the algorithm, determines at
         what point an image should not be updated further.
         Default is 20% of the FDK L2 norm.
--------------------------------------------------------------------
    This file is part of the TIGRE Toolbox

    Copyright (c) 2015, University of Bath and
                        CERN-European Organization for Nuclear Research
                        All rights reserved.

    License:            Open Source under BSD.
                        See the full license at
                        https://github.com/CERN/TIGRE/license.txt

    Contact:            tigre.toolbox@gmail.com
    Codes:              https://github.com/CERN/TIGRE/
    --------------------------------------------------------------------
    Coded by:          MATLAB (original code): Ander Biguri
                       PYTHON : Reuben Lindroos



    """

    def __init__(self, proj, geo, angles, niter, **kwargs):
        
        # if "blocksize" not in kwargs:
        #     kwargs.update(dict(blocksize=1))
        #kwargs.update(dict(regularisation="minimizeTV"))
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        if "maxl2err" not in kwargs:
            self.epsilon = (
                im3DNORM(Ax(FDK(proj, geo, angles, gpuids=self.gpuids), geo, angles) - proj, 2)
                * 0.2
            )
        else:
            self.epsilon = kwargs["maxl2err"]
        self.numiter_tv = 20 if "tviter" not in kwargs else kwargs["tviter"]
        self.beta = self.lmbda
        self.beta_red = self.lmbda_red

    # Override
    def run_main_iter(self):
        stop_criteria = False
        n_iter = 0
        d_p_1st = 1
        while not stop_criteria:
            if self.verbose:
                self._estimate_time_until_completion(n_iter)
                
            res_prev = copy.deepcopy(self.res)
            n_iter += 1
            est_proj = Ax(self.res, self.geo, self.angles)
            d_p = im3DNORM(est_proj - self.proj, 2)
            if d_p**2 > self.epsilon:
                getattr(self, self.dataminimizing)()
            dd = im3DNORM(Ax(self.res,self.geo,self.angles)-self.proj, 2)
            dp_vec = self.res - res_prev
            dp = im3DNORM(dp_vec, 2)
            
            dtvg = 1 if n_iter == 1 else d_p / d_p_1st
                
            res_prev = copy.deepcopy(self.res)
            self.res = getattr(self, self.regularisation)(self.res, dtvg)
            dg_vec = self.res - res_prev
            dg = im3DNORM(dg_vec, 2)
            
            if n_iter == 1:
                d_p_1st = im3DNORM(Ax(res_prev,self.geo,self.angles)-self.proj, 2)

            self.beta *= self.beta_red
            c = np.dot(dg_vec.reshape(-1,), dp_vec.reshape(-1,)) / max(dg * dp, 1e-6) # reshape ensures no copy is made. 
            if (c < -0.99 and dd <=
                    self.epsilon) or self.beta < 0.005 or n_iter > self.niter:
                if self.verbose:
                    print("\n"
                          "     Stop criteria met: \n"
                          "     c = " + str(c) + "\n"
                          "     beta = " + str(self.beta) + "\n"
                          "     iter = " + str(n_iter) + "\n")
                stop_criteria = True


class PCSD(initPCSD):
    __doc__ = initPCSD.__doc__
    
    def __init__(self, proj, geo, angles, niter, **kwargs):

        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to 1, please use an OS version of the algorithm for blocksize > 1')
        kwargs.update(dict(blocksize=1))
        kwargs.update(dict(regularisation="minimizeTV"))
        
        initPCSD.__init__(self, proj, geo, angles, niter, **kwargs)
        
    
pcsd = decorator(PCSD, name='pcsd')


class AwPCSD(initPCSD): 
    __doc__ = (        
        "    Adaptive Weighted TV (edge preserving) version of PCSD\n\n"
        "    :extra keyword delta: (float)\n"
        "        Controls amount of smoothing at edges of the image \n"
        "        Default -0.005\n"
        ) + initPCSD.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        if "blocksize" in kwargs and kwargs['blocksize']>1:
            print('Warning: blocksize is set to 1, please use an OS version of the algorithm for blocksize > 1')
        kwargs.update(dict(blocksize=1))
        kwargs.update(dict(regularisation="minimizeAwTV"))
        self.delta = np.float32(-0.005) if "delta" not in kwargs else kwargs["delta"]
            
        initPCSD.__init__(self, proj, geo, angles, niter, **kwargs)
        

aw_pcsd = decorator(AwPCSD, name="aw_pcsd")


class OS_PCSD(initPCSD):
    __doc__ = (
        "    Oriented Subsets version of PCSD\n\n"
        "        Default blocksize = 20"
        ) + initPCSD.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        self.blocksize = 20 if "blocksize" not in kwargs else kwargs["blocksize"]
        kwargs.update(dict(regularisation="minimizeTV"))
            
        initPCSD.__init__(self, proj, geo, angles, niter, **kwargs)
        

os_pcsd = decorator(OS_PCSD, name="os_pcsd")


class OS_Aw_PCSD(initPCSD): 
    __doc__ = (
        "    Oriented Subsets and Adaptive Weighted TV version of PCSD\n\n"
        "    Default blocksize = 20\n\n"
        "    :extra keyword delta: (float)\n"
        "        Controls amount of smoothing at edges of the image \n"
        "        Default -0.005\n"
        ) + initPCSD.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        self.blocksize = 20 if "blocksize" not in kwargs else kwargs["blocksize"]
        kwargs.update(dict(regularisation="minimizeAwTV"))
        self.delta = np.float32(-0.005) if "delta" not in kwargs else kwargs["delta"]
            
        initPCSD.__init__(self, proj, geo, angles, niter, **kwargs)
        
        
os_aw_pcsd = decorator(OS_Aw_PCSD, name="os_aw_pcsd")



