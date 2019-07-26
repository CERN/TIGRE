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
        Reduction of lambda every iteration
        lambda=lambdared*lambda. Default is 0.99

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

        if 'blocksize' not in kwargs:
            kwargs.update(blocksize=1)
        IterativeReconAlg.__init__(self, proj, geo, angles, niter, **kwargs)
        if 'alpha' not in kwargs:
            self.alpha = 0.002
        if 'alpha_red' not in kwargs:
            self.alpha_red = 0.95
        if 'rmax' not in kwargs:
            self.rmax = 0.95
        if 'maxl2err' not in kwargs:
            self.epsilon = im3DNORM(FDK(proj, geo, angles), 2) * 0.2
        else:
            self.epsilon = kwargs['maxl2err']
        if "tviter" not in kwargs:
            self.numiter_tv = 20
        else:
            self.numiter_tv = kwargs["tviter"]
        if 'regularisation' not in kwargs:
            self.regularisation = 'minimizeTV'
        self.beta = self.lmbda
        self.beta_red = self.lmbda_red

    # Override
    def run_main_iter(self):
        stop_criteria = False
        n_iter = 0
        while not stop_criteria:
            if self.verbose:
                if n_iter == 0:
                    print("POCS Algorithm in progress.")
                    toc = time.clock()
                if n_iter == 1:
                    tic = time.clock()
                    print('Esitmated time until completetion (s): ' +
                          str((self.niter - 1) * (tic - toc)))
            res_prev = copy.deepcopy(self.res)
            n_iter += 1
            getattr(self, self.dataminimizing)()
            g = Ax(self.res, self.geo, self.angles)
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
            c = np.dot(dg_vec.reshape(-1,), dp_vec.reshape(-1,)) / \
                max(dg * dp, 1e-6)
            if (c < -0.99 and dd <=
                    self.epsilon) or self.beta < 0.005 or n_iter > self.niter:
                if self.verbose:
                    print("\n"
                          "     Stop criteria met: \n"
                          "     c = " + str(c) + "\n"
                          "     beta = " + str(self.beta) + "\n"
                          "     iter = " + str(n_iter) + "\n")
                stop_criteria = True


asd_pocs = decorator(ASD_POCS, name='asd_pocs')


class AwASD_POCS(ASD_POCS):
    __doc__ = ASD_POCS.__doc__

    def __init__(self, proj, geo, angles, niter, **kwargs):

        kwargs.update(dict(regularisation='minimizeAwTV'))

        if 'delta' not in kwargs:
            self.delta = np.array([-0.005], dtype=np.float32)[0]
        ASD_POCS.__init__(self, proj, geo, angles, niter, **kwargs)


awasd_pocs = decorator(AwASD_POCS, name='awasd_pocs')
