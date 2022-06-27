from __future__ import division

import copy
import time

import numpy as np
from _AwminTV import AwminTV
from _minTV import minTV
from tigre.algorithms.single_pass_algorithms import FDK
from tigre.utilities.Atb import Atb
from tigre.utilities.Ax import Ax
from tigre.utilities.im3Dnorm import im3DNORM
from tigre.utilities.init_multigrid import init_multigrid
from tigre.utilities.order_subsets import order_subsets
from tigre.utilities.Measure_Quality import Measure_Quality as MQ
from tigre.utilities.gpu import GpuIds

"""
This module is where the umbrella class IterativeReconAlg is located
which is the umbrella class to all the other algorithms apart from
the single pass type algorithms.
"""

# coding: utf8


if hasattr(time, "perf_counter"):
    default_timer = time.perf_counter
else:
    default_timer = time.clock


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

        self.geo.check_geo(angles)

        options = dict(
            blocksize=20,
            lmbda=1,
            lmbda_red=1,
            OrderStrategy=None,
            Quameasopts=None,
            init=None,
            verbose=True,
            noneg=True,
            computel2=False,
            dataminimizing="art_data_minimizing",
            name="Iterative Reconstruction",
            sup_kw_warning=False,
            gpuids=None,
        )
        allowed_keywords = [
            "V",
            "W",
            "log_parameters",
            "angleblocks",
            "angle_index",
            "alpha",
            "alpha_red",
            "rmax",
            "maxl2err",
            "delta",
            "regularisation",
            "tviter",
            "tvlambda",
            "hyper",
            "fista_p",
            "fista_q"
        ]
        self.__dict__.update(options)
        self.__dict__.update(**kwargs)
        for kw in kwargs.keys():
            if kw not in options and (kw not in allowed_keywords):
                if self.verbose:
                    if not kwargs.get("sup_kw_warning"):
                        # Note: might not want this warning (typo checking).
                        print(
                            "Warning: "
                            + kw
                            + " not recognised as default parameter for instance of IterativeReconAlg."  # noqa: E501
                        )
        if self.angles.ndim == 1:
            a1 = self.angles
            a2 = np.zeros(self.angles.shape[0], dtype=np.float32)
            self.angles = np.vstack((a1, a2, a2)).T
        if not all([hasattr(self, "angleindex"), hasattr(self, "angleblocks")]):
            self.set_angle_index()
        if not hasattr(self, "W"):
            self.set_w()
        if not hasattr(self, "V"):
            self.set_v()
        if not hasattr(self, "res"):
            self.set_res()
        if self.verbose:
            self.tic = 0    # preparation for _estimate_time_until_completion()
        # make it list
        if self.Quameasopts is not None:
            self.Quameasopts = (
                [self.Quameasopts] if isinstance(self.Quameasopts, str) else self.Quameasopts
            )
            setattr(self, "lq", np.zeros([len(self.Quameasopts), niter]))  # quameasoptslist
        else:
            setattr(self, "lq", np.zeros([0, niter]))  # quameasoptslist
        setattr(self, "l2l", np.zeros([1, niter]))  # l2list

    def set_w(self):
        """
        Calculates value of W if this is not given.
        :return: None
        """
        geox = copy.deepcopy(self.geo)
        geox.sVoxel[1:] = geox.sVoxel[1:] * 1.1  # a bit larger to avoid zeros in projections
        geox.sVoxel[0] = max(geox.sDetector[0], geox.sVoxel[0])

        geox.nVoxel = np.array([2, 2, 2])
        geox.dVoxel = geox.sVoxel / geox.nVoxel
        W = Ax(
            np.ones(geox.nVoxel, dtype=np.float32), geox, self.angles, "Siddon", gpuids=self.gpuids
        )
        W[W <= min(self.geo.dVoxel / 2)] = np.inf
        W = 1.0 / W
        setattr(self, "W", W)

    def set_v(self):
        """
        Computes value of V parameter if this is not given.
        :return: None
        """
        block_count = len(self.angleblocks)
        geo = self.geo
        V = np.ones((block_count, geo.nVoxel[1], geo.nVoxel[2]), dtype=np.float32)

        for i in range(block_count):
            if geo.mode != "parallel":

                geox = copy.deepcopy(self.geo)
                geox.angles = self.angleblocks[i]

                geox.DSD = geo.DSD[self.angle_index[i]]
                geox.DSO = geo.DSO[self.angle_index[i]]
                geox.offOrigin = geo.offOrigin[self.angle_index[i], :]
                geox.offDetector = geo.offDetector[self.angle_index[i], :]
                geox.rotDetector = geo.rotDetector[self.angle_index[i], :]
                geox.COR = geo.COR[self.angle_index[i]]

                # shrink the volume size to avoid zeros in backprojection
                geox.sVoxel = (
                    geox.sVoxel * np.max(geox.sVoxel[1:] / np.linalg.norm(geox.sVoxel[1:])) * 0.9
                )
                geox.dVoxel = geox.sVoxel / geox.nVoxel
                proj_one = np.ones(
                    (len(self.angleblocks[i]), geo.nDetector[0], geo.nDetector[1]), dtype=np.float32
                )
                V[i] = Atb(proj_one, geox, self.angleblocks[i], "FDK", gpuids=self.gpuids).mean(
                    axis=0
                )

            else:
                V[i] *= len(self.angleblocks[i])
        V[V==0.0] = np.inf       

        self.V = V

    def set_res(self):
        """
        Calulates initial value for res if this is not given.
        :return: None
        """
        self.res = np.zeros(self.geo.nVoxel, dtype=np.float32)
        init = self.init
        verbose = self.verbose
        if isinstance(init, str):
            if init == "multigrid":
                if verbose:
                    print("init multigrid in progress...")
                    print("default blocksize=1 for init_multigrid(OS_SART)")
                self.res = init_multigrid(self.proj, self.geo, self.angles, alg="SART")
                if verbose:
                    print("init multigrid complete.")
            if init == "FDK":
                self.res = FDK(self.proj, self.geo, self.angles)

        elif isinstance(init, np.ndarray):
            if (self.geo.nVoxel == init.shape).all():
                self.res = init
            else:
                raise ValueError("wrong dimension of array for initialisation")
        elif init is not None:
            raise ValueError("wrong value for initialisation")

    def set_angle_index(self):
        """
        sets angle_index and angleblock if this is not given.
        :return: None
        """
        self.angleblocks, self.angle_index = order_subsets(
            self.angles, self.blocksize, self.OrderStrategy
        )

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
                self._estimate_time_until_completion(i)

            getattr(self, self.dataminimizing)()
            self.error_measurement(res_prev, i)

    def art_data_minimizing(self):
        geo = copy.deepcopy(self.geo)

        for j in range(len(self.angleblocks)):

            if self.blocksize == 1:
                angle = np.array([self.angleblocks[j]], dtype=np.float32)
                angle_indices = np.array([self.angle_index[j]], dtype=np.int32)

            else:
                angle = self.angleblocks[j]
                angle_indices = self.angle_index[j]
                # slice parameters if needed
            geo.offOrigin = self.geo.offOrigin[angle_indices]
            geo.offDetector = self.geo.offDetector[angle_indices]
            geo.rotDetector = self.geo.rotDetector[angle_indices]
            geo.DSD = self.geo.DSD[angle_indices]
            geo.DSO = self.geo.DSO[angle_indices]

            self.update_image(geo, angle, j)

            if self.noneg:
                self.res = self.res.clip(min=0)

    def minimizeTV(self, res_prev, dtvg):
        if self.gpuids is None:
            self.gpuids = GpuIds()
        return minTV(res_prev, dtvg, self.numiter_tv, self.gpuids)

    def minimizeAwTV(self, res_prev, dtvg):
        if self.gpuids is None:
            self.gpuids = GpuIds()
        return AwminTV(res_prev, dtvg, self.numiter_tv, self.delta, self.gpuids)

    def error_measurement(self, res_prev, iter):
        if self.Quameasopts is not None:
            self.lq[:, iter] = MQ(self.res, res_prev, self.Quameasopts)
        if self.computel2:
            # compute l2 borm for b-Ax
            errornow = im3DNORM(
                self.proj - Ax(self.res, self.geo, self.angles, "Siddon", gpuids=self.gpuids), 2
            )
            self.l2l[0, iter] = errornow

    def update_image(self, geo, angle, iteration):
        """
        VERBOSE:
         for j in range(angleblocks):
             angle = np.array([alpha[j]], dtype=np.float32)
             proj_err = proj[angle_index[j]] - Ax(res, geo, angle, 'Siddon')
             weighted_err = W[angle_index[j]] * proj_err
             backprj = Atb(weighted_err, geo, angle, 'FDK')
             weighted_backprj = 1 / V[angle_index[j]] * backprj
             res += weighted_backprj
             res[res<0]=0

        :return: None
        """

        ang_index = self.angle_index[iteration].astype(np.int)

        self.res += (
            self.lmbda
            * 1.0
            / self.V[iteration]
            * Atb(
                self.W[ang_index]
                * (self.proj[ang_index] - Ax(self.res, geo, angle, "Siddon", gpuids=self.gpuids)),
                geo,
                angle,
                "FDK",
                gpuids=self.gpuids,
            )
        )

    def getres(self):
        return self.res

    def geterrors(self):
        if self.computel2:
            return np.concatenate((self.l2l, self.lq), axis=0)
        else:
            return self.lq

    def __str__(self):
        parameters = []
        for item in self.__dict__:
            if item == "geo":
                pass
            elif hasattr(self.__dict__.get(item), "shape"):
                if self.__dict__.get(item).ravel().shape[0] > 100:
                    parameters.append(item + " shape: " + str(self.__dict__.get(item).shape))
            else:
                parameters.append(item + ": " + str(self.__dict__.get(item)))

        return "\n".join(parameters)

    def _estimate_time_until_completion(self, iter):
        if iter == 0:
            print(str(self.name).upper() + " " + "algorithm in progress.")
            self.tic = default_timer()
        if iter == 1:
            toc = default_timer()

            remaining_time = (self.niter - 1) * (toc - self.tic)
            seconds = int(remaining_time)
            print(
                "Estimated time until completion : "
                + time.strftime("%H:%M:%S", time.gmtime(seconds))
            )


def decorator(IterativeReconAlg, name=None, docstring=None):  # noqa: N803
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
    >>> geo = tigre.geometry_defaut(high_resolution=False)
    >>> src = load_head_phantom(number_of_voxels=geo.nVoxel)
    >>> proj = Ax(src,geo,angles)
    >>> angles = np.linspace(0,2*np.pi,100)
    >>> iterativereconalg = decorator(IterativeReconAlg)
    >>> output = iterativereconalg(proj,geo,angles, niter=50)

    """

    def iterativereconalg(proj, geo, angles, niter, **kwargs):

        geo.check_geo(angles)
        alg = IterativeReconAlg(proj, geo, angles, niter, **kwargs)
        if name is not None:
            alg.name = name
        alg.run_main_iter()
        if alg.computel2 or alg.Quameasopts is not None:
            return alg.getres(), alg.geterrors()
        else:
            return alg.getres()

    if docstring is not None:
        iterativereconalg.__doc__ = docstring + IterativeReconAlg.__doc__
    else:
        iterativereconalg.__doc__ = IterativeReconAlg.__doc__
    if name is not None:
        iterativereconalg.__name__ = name

    return iterativereconalg


iterativereconalg = decorator(IterativeReconAlg)
