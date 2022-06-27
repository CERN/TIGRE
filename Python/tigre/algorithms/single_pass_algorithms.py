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

    proj_filt = filtering(proj_filt, geo, angles, parker=False, verbose=verbose)
    
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
    proj_filt = filtering(copy.deepcopy(proj), geox, angles, parker=False, verbose=verbose)
    return Atb(proj_filt, geo, angles, gpuids=gpuids) * geo.DSO / geo.DSD
