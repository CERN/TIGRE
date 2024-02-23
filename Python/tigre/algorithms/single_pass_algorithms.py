from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy

import numpy as np
from tigre.utilities.Atb import Atb
from tigre.utilities.filtering import filtering


def FDK(proj, geo, angles, **kwargs):
    """
    solves CT image reconstruction.

    :param proj: np.array(dtype=float32),
    Data input in the form of 3d

    :param geo: tigre.utilities.geometry.Geometry
    Geometry of detector and image (see examples/Demo code)

    :param angles: np.array(dtype=float32)
    Angles of projection, shape = (nangles,3) or (nangles,)

    :param filter: str
    Type of filter used for backprojection
    opts: "shep_logan"
          "cosine"
          "hamming"
          "hann"

    :param verbose: bool
    Feedback print statements for algorithm progress

    :param kwargs: dict
    keyword arguments

    :return: np.array(dtype=float32)

    Usage:
    -------
    >>> import tigre
    >>> import tigre.algorithms as algs
    >>> import numpy
    >>> from tigre.demos.Test_data import data_loader
    >>> geo = tigre.geometry(mode='cone',default_geo=True,
    >>>                         nVoxel=np.array([64,64,64]))
    >>> angles = np.linspace(0,2*np.pi,100)
    >>> src_img = data_loader.load_head_phantom(geo.nVoxel)
    >>> proj = tigre.Ax(src_img,geo,angles)
    >>> output = algs.FDK(proj,geo,angles)

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
    Coded by:           MATLAB (original code): Ander Biguri
                        PYTHON : Reuben Lindroos
    """
    verbose = kwargs["verbose"] if "verbose" in kwargs else False

    gpuids = kwargs["gpuids"] if "gpuids" in kwargs else None
    dowang = kwargs["dowang"] if "dowang" in kwargs else False

    def zeropadding(proj, geo):
        zgeo = copy.deepcopy(geo)
        padwidth = int(2 * geo.offDetector[1] / geo.dDetector[1])
        zgeo.offDetector[1] = geo.offDetector[1] - \
            padwidth / 2 * geo.dDetector[1]
        zgeo.nDetector[1] = abs(padwidth) + geo.nDetector[1]
        zgeo.sDetector[1] = zgeo.nDetector[1] * zgeo.dDetector[1]

        theta = (geo.sDetector[1] / 2 - abs(geo.offDetector[1])
                 ) * np.sign(geo.offDetector[1])

        print(padwidth, zgeo.nDetector[1], zgeo.offDetector[1], theta)
        if geo.offDetector[1] > 0:
            zproj = np.zeros(
                (proj.shape[0] , proj.shape[1], proj.shape[2]+ padwidth), dtype=proj.dtype)
            for ii in range(proj.shape[0]):
                zproj[ii,:, :] = np.concatenate(
                    (np.zeros((proj.shape[1], padwidth)), proj[:, :, ii]), axis=1)
        else:
            zproj = np.zeros(
                (proj.shape[0] , proj.shape[1] , proj.shape[2]+ abs(padwidth)), dtype=proj.dtype)
            for ii in range(proj.shape[0]):
                zproj[ii, :, :] = np.concatenate(
                    (proj[ ii,:, :], np.zeros((proj.shape[1], abs(padwidth)))), axis=1)

        return zproj, zgeo, theta

    def preweighting2(proj, geo, theta):
        """
        Preweighting using Wang function

        :param proj: np.array(dtype=float32),
        Data input in the form of 3d

        :param geo: tigre.utilities.geometry.Geometry
        Geometry of detector and image (see examples/Demo code)

        :param theta: np.array(dtype=float32)
        Angles of projection, shape = (nangles,3) or (nangles,)

        :return: np.array(dtype=float32), np.array(dtype=float32)

        """
        offset = geo.offDetector[1]
        offset = offset + (geo.DSD / geo.DSO) #* geo.COR

        us = np.arange(-geo.nDetector[1]/2 + 0.5, geo.nDetector[1] /
                       2 - 0.5 + 1) * geo.dDetector[1] + abs(offset)
        us = us * geo.DSO / geo.DSD
        abstheta = np.abs(theta * geo.DSO / geo.DSD)

        w = np.ones(proj[0, :, :].shape)

        for ii in range(geo.nDetector[1]):
            t = us[ii]
            if np.abs(t) <= abstheta:
                w[:,ii] = 0.5 * (np.sin((np.pi / 2) * np.arctan(t /
                                  geo.DSO) / (np.arctan(abstheta / geo.DSO))) + 1)
            if t < -abstheta:
                w[:,ii] = 0

        if theta < 0:
            w = np.fliplr(w)

        proj_w = proj.copy()  # preallocation
        for ii in range(proj.shape[0]):
            proj_w[ii, :, :,] = proj[ii, :, :] * w * 2

        return proj_w, w


    if dowang:
        print('FDK: applying detector offset weights')
        # Zero-padding to avoid FFT-induced aliasing
        zproj, zgeo, theta = zeropadding(proj, geo)
        # Preweighting using Wang function to save memory
        proj, _ = preweighting2(zproj, zgeo, theta)

        # Replace original proj and geo
        # proj = proj_w;
        geo = zgeo

    print(geo)
    print(proj.shape)
    geo = copy.deepcopy(geo)
    geo.check_geo(angles)
    geo.checknans()
    geo.filter = kwargs["filter"] if "filter" in kwargs else None
    # Weight
    proj_filt = np.zeros(proj.shape, dtype=np.float32)
    xv = np.arange((-geo.nDetector[1] / 2) + 0.5,
                   1 + (geo.nDetector[1] / 2) - 0.5) * geo.dDetector[1]
    yv = np.arange((-geo.nDetector[0] / 2) + 0.5,
                   1 + (geo.nDetector[0] / 2) - 0.5) * geo.dDetector[0]
    (yy, xx) = np.meshgrid(xv, yv)

    w = geo.DSD[0] / np.sqrt((geo.DSD[0] ** 2 + xx ** 2 + yy ** 2))
    np.multiply(proj, w, out=proj_filt)

    proj_filt = filtering(proj_filt, geo, angles,
                          parker=False, verbose=verbose)
    # geo.proj = proj_filt

    return Atb(proj_filt, geo, geo.angles, "FDK", gpuids=gpuids)


fdk = FDK


def fbp(proj, geo, angles, **kwargs):  # noqa: D103
    __doc__ = FDK.__doc__  # noqa: F841
    if geo.mode != "parallel":
        raise ValueError("Only use FBP for parallel beam. Check geo.mode.")
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    verbose = kwargs["verbose"] if "verbose" in kwargs else False
    gpuids = kwargs["gpuids"] if "gpuids" in kwargs else None
    proj_filt = filtering(copy.deepcopy(proj), geox,
                          angles, parker=False, verbose=verbose)
    return Atb(proj_filt, geo, angles, gpuids=gpuids) * geo.DSO / geo.DSD
