from __future__ import division

import numpy as np
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
from tigre.utilities.order_subsets import order_subsets
from tigre.utilities.init_multigrid import init_multigrid
from tigre.utilities.Measure_Quality import Measure_Quality as MQ
from tigre.utilities.im3Dnorm import im3DNORM
from tigre.algorithms.single_pass_algorithms import FDK
from _minTV import minTV
from _AwminTV import AwminTV
import time
import copy
"""
This module is where the umbrella class IterativeReconAlg is located
which is the umbrella class to all the other algorithms apart from
the single pass type algorithms.
"""

# coding: utf8


class IterativeReconAlg(object):
    """
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
        of proj and backproj for OS_SART
    :keyword lmbda: (np.float64)
        Sets the value of the hyperparameter.

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
        For algorithms that make use of a tvdenoising step in their
        iterations. This includes:

            OS_SART_TV
            ASD_POCS
            AWASD_POCS
            FISTA

    :keyword tvlambda: (float)
        For algorithms that make use of a tvdenoising step in their
        iterations.

            OS_SART_TV
            FISTA

    Usage
    --------
    >>> import numpy as np
    >>> import tigre
    >>> import tigre.algorithms as algs
    >>> from tigre.demos.Test_data import data_loader
    >>> geo = tigre.geometry(mode='cone',default_geo=True,
    >>>                         nVoxel=np.array([64,64,64]))
    >>> angles = np.linspace(0,2*np.pi,100)
    >>> src_img = data_loader.load_head_phantom(geo.nVoxel)
    >>> proj = tigre.Ax(src_img,geo,angles)
    >>> output = algs.iterativereconalg(proj,geo,angles,niter=50
    >>>                                 blocksize=20)

    tigre.demos.run() to launch ipython notebook file with examples.


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

        self.proj = proj
        self.angles = angles
        self.geo = geo
        self.niter = niter

        options = dict(blocksize=20, lmbda=1, lmbda_red=0.99,
                       OrderStrategy=None, Quameasopts=None,
                       init=None, verbose=True, noneg=True,
                       computel2=False, dataminimizing='art_data_minimizing',
                       name='Iterative Reconstruction', sup_kw_warning=False)
        allowed_keywords = [
            'V',
            'W',
            'log_parameters',
            'angleblocks',
            'angle_index',
            'delta',
            'regularisation',
            'tviter',
            'tvlambda',
            'hyper']
        self.__dict__.update(options)
        self.__dict__.update(**kwargs)
        for kw in kwargs.keys():
            if kw not in options and (kw not in allowed_keywords):
                if self.verbose:
                    if not kwargs.get('sup_kw_warning'):
                        # Note: might not want this warning (typo checking).
                        print(
                            "Warning: " +
                            kw +
                            " not recognised as default parameter for instance of IterativeReconAlg.")
        if self.angles.ndim == 1:
            a1 = self.angles
            a2 = np.zeros(self.angles.shape[0], dtype=np.float32)
            setattr(self, 'angles', np.vstack((a1, a2, a2)).T)
        if not all([hasattr(self, 'angleindex'),
                    hasattr(self, 'angleblocks')]):
            self.set_angle_index()
        if not hasattr(self, 'W'):
            self.set_w()
        if not hasattr(self, 'V'):
            self.set_v()
        if not hasattr(self, 'res'):
            self.set_res()
        setattr(self, 'lq', [])  # quameasoptslist
        setattr(self, 'l2l', [])  # l2list

    def set_w(self):
        """
        Calculates value of W if this is not given.
        :return: None
        """
        geox = copy.deepcopy(self.geo)
        geox.sVoxel[0:] = self.geo.DSD - self.geo.DSO
        geox.sVoxel[2] = max(geox.sDetector[1], geox.sVoxel[2])
        geox.nVoxel = np.array([2, 2, 2])
        geox.dVoxel = geox.sVoxel / geox.nVoxel
        W = Ax(
            np.ones(
                geox.nVoxel,
                dtype=np.float32),
            geox,
            self.angles,
            "ray-voxel")
        W[W <= min(self.geo.dVoxel / 4)] = np.inf
        W = 1. / W
        setattr(self, 'W', W)

    def set_v(self):
        """
        Computes value of V parameter if this is not given.
        :return: None
        """
        geo = self.geo
        if geo.mode != 'parallel':

            start = geo.sVoxel[1] / 2 - geo.dVoxel[1] / 2 + geo.offOrigin[1]
            stop = -geo.sVoxel[1] / 2 + geo.dVoxel[1] / 2 + geo.offOrigin[1]
            step = -geo.dVoxel[1]

            xv = np.arange(start, stop + step, step)

            start = geo.sVoxel[2] / 2 - geo.dVoxel[2] / 2 + geo.offOrigin[2]
            stop = -geo.sVoxel[2] / 2 + geo.dVoxel[2] / 2 + geo.offOrigin[2]
            step = -geo.dVoxel[2]

            yv = -1 * np.arange(start, stop + step, step)

            (yy, xx) = np.meshgrid(yv, xv)
            A = (self.angles[:, 0] + np.pi / 2)

            V = np.empty((self.angles.shape[0], geo.nVoxel[1], geo.nVoxel[2]))
            for i in range(self.angles.shape[0]):
                if hasattr(geo.DSO, 'shape') and len(geo.DSO.shape) >= 1:
                    DSO = geo.DSO[i]
                else:
                    DSO = geo.DSO
                V[i] = (DSO / (DSO + (yy * np.sin(-A[i])) -
                               (xx * np.cos(-A[i])))) ** 2

        else:
            V = np.ones((self.angles.shape[0]), dtype=np.float32)
        if self.blocksize > 1:
            v_list = [np.sum(V[self.angle_index[i]], axis=0)
                      for i in range(len(self.angleblocks))]
            V = np.stack(v_list, 0)

        V = np.array(V, dtype=np.float32)
        setattr(self, 'V', V)

    def set_res(self):
        """
        Calulates initial value for res if this is not given.
        :return: None
        """
        setattr(self, 'res', np.zeros(self.geo.nVoxel, dtype=np.float32))
        init = self.init
        verbose = self.verbose
        if init == 'multigrid':
            if verbose:
                print('init multigrid in progress...')
                print('default blocksize=1 for init_multigrid(OS_SART)')
            self.res = init_multigrid(
                self.proj, self.geo, self.angles, alg='SART')
            if verbose:
                print('init multigrid complete.')
        if init == 'FDK':
            self.res = FDK(self.proj, self.geo, self.angles)

        if isinstance(init, np.ndarray):
            if (self.geo.nVoxel == init.shape).all():

                self.res = init

            else:
                raise ValueError('wrong dimension of array for initialisation')

    def set_angle_index(self):
        """
        sets angle_index and angleblock if this is not given.
        :return: None
        """
        angleblocks, angle_index = order_subsets(
            self.angles, self.blocksize, self.OrderStrategy)
        setattr(self, 'angleblocks', angleblocks)
        setattr(self, 'angle_index', angle_index)

    def run_main_iter(self):
        """
        Goes through the main iteration for the given configuration.
        :return: None
        """
        Quameasopts = self.Quameasopts

        for i in range(self.niter):

            res_prev = None
            if Quameasopts is not None:
                res_prev = copy.deepcopy(self.res)
            if self.verbose:
                if i == 0:
                    print(str(self.name).upper() +
                          ' ' + "algorithm in progress.")
                    toc = time.clock()
                if i == 1:
                    tic = time.clock()
                    print('Esitmated time until completetion (s): ' +
                          str((self.niter - 1) * (tic - toc)))
            getattr(self, self.dataminimizing)()
            self.error_measurement(res_prev, i)

    def art_data_minimizing(self):

        geo = copy.deepcopy(self.geo)
        for j in range(len(self.angleblocks)):
            if self.blocksize == 1:
                angle = np.array([self.angleblocks[j]], dtype=np.float32)
            else:
                angle = self.angleblocks[j]

            if geo.offOrigin.shape[0] == self.angles.shape[0]:
                geo.offOrigin = self.geo.offOrigin[j]
            if geo.offDetector.shape[0] == self.angles.shape[0]:
                geo.offOrin = self.geo.offDetector[j]
            if geo.rotDetector.shape[0] == self.angles.shape[0]:
                geo.rotDetector = self.geo.rotDetector[j]
            if hasattr(geo.DSD, 'shape') and len((geo.DSD.shape)):
                if geo.DSD.shape[0] == self.angles.shape[0]:
                    geo.DSD = self.geo.DSD[j]
            if hasattr(geo.DSO, 'shape') and len((geo.DSD.shape)):
                if geo.DSO.shape[0] == self.angles.shape[0]:
                    geo.DSO = self.geo.DSO[j]

            self.update_image(geo, angle, j)

            if self.noneg:
                self.res = self.res.clip(min=0)


    def minimizeTV(self, res_prev, dtvg):
        return minTV(res_prev, dtvg, self.numiter_tv)

    def minimizeAwTV(self, res_prev, dtvg):
        return AwminTV(res_prev, dtvg, self.numiter_tv, self.delta)

    def error_measurement(self, res_prev, iter):
        if self.Quameasopts is not None and iter > 0:
            self.lq.append(MQ(self.res, res_prev, self.Quameasopts))
        if self.computel2:
            # compute l2 borm for b-Ax
            errornow = im3DNORM(
                self.proj - Ax(self.res, self.geo, self.angles, 'ray-voxel'), 2)
            self.l2l.append(errornow)

    def update_image(self, geo, angle, iteration):
        """
        VERBOSE:
         for j in range(angleblocks):
             angle = np.array([alpha[j]], dtype=np.float32)
             proj_err = proj[angle_index[j]] - Ax(res, geo, angle, 'ray-voxel')
             weighted_err = W[angle_index[j]] * proj_err
             backprj = Atb(weighted_err, geo, angle, 'FDK')
             weighted_backprj = 1 / V[angle_index[j]] * backprj
             res += weighted_backprj
             res[res<0]=0

        :return: None
        """
        self.res += self.lmbda * 1. / self.V[iteration] * Atb(self.W[self.angle_index[iteration]] * (
            self.proj[self.angle_index[iteration]] - Ax(self.res, geo, angle, 'interpolated')), geo, angle, 'FDK')

    def getres(self):
        return self.res

    def geterrors(self):
        return self.l2l, self.lq

    def __str__(self):
        parameters = []
        for item in self.__dict__:
            if item == 'geo':
                pass
            elif hasattr(self.__dict__.get(item), 'shape'):
                if self.__dict__.get(item).ravel().shape[0] > 100:
                    parameters.append(item + ' shape: ' +
                                      str(self.__dict__.get(item).shape))
            else:
                parameters.append(item + ': ' + str(self.__dict__.get(item)))

        return '\n'.join(parameters)


def decorator(IterativeReconAlg, name=None, docstring=None):
    """
    Calls run_main_iter when parameters are given to it.

    :param IterativeReconAlg: obj, class
        instance of IterativeReconAlg
    :param name: str
        for name of func
    :param docstring: str
        other documentation that may need to be included from external source.
    :return: func

    Examples
    --------
    >>> import tigre
    >>> from tigre.demos.Test_data.data_loader import load_head_phantom
    >>> geo = tigre.geometry_defaut(high_quality=False)
    >>> src = load_head_phantom(number_of_voxels=geo.nVoxel)
    >>> proj = Ax(src,geo,angles)
    >>> angles = np.linspace(0,2*np.pi,100)
    >>> iterativereconalg = decorator(IterativeReconAlg)
    >>> output = iterativereconalg(proj,geo,angles, niter=50)

    """

    def iterativereconalg(proj, geo, angles, niter, **kwargs):
        alg = IterativeReconAlg(proj, geo, angles, niter, **kwargs)
        if name is not None:
            alg.name = name
        alg.run_main_iter()
        if alg.computel2:
            return alg.getres(), alg.geterrors()
        else:
            return alg.getres()

    if docstring is not None:
        setattr(
            iterativereconalg,
            '__doc__',
            docstring +
            IterativeReconAlg.__doc__)
    else:
        setattr(iterativereconalg, '__doc__', IterativeReconAlg.__doc__)
    if name is not None:
        setattr(iterativereconalg, '__name__', name)
    return iterativereconalg


iterativereconalg = decorator(IterativeReconAlg)
