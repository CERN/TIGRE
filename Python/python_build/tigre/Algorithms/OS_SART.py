from __future__ import division
from __future__ import print_function
import time
import os
import sys
from tigre.Utilities.init_multigrid import init_multigrid
from scipy.linalg import *
import numpy as np
import copy
from _Ax import Ax
from _Atb import Atb
from tigre.Utilities.order_subsets import order_subsets
from tigre.Utilities.Measure_Quality import Measure_Quality as MQ
from tigre.Algorithms.FDK import FDK
from tigre.Utilities.im3Dnorm import im3DNORM
# TODO: this is quite nasty; it would be nice to reorganise file structure later so top level folder is always in path
currDir = os.path.dirname(os.path.realpath(__file__))
rootDir = os.path.abspath(os.path.join(currDir, '..'))
if rootDir not in sys.path:  # add parent dir to paths
    sys.path.append(rootDir)

from scipy.linalg import *


def OS_SART(proj, geo, alpha, niter,
            blocksize=20, lmbda=1, lmbda_red=0.99, OrderStrategy=None, Quameasopts=None, init=None, verbose=True,noneg=True,computel2=False):
    ('\n'
     """SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets
              Simultaneous Algebraic Reconxtruction Techique algorithm

   OS_SART(PROJ,GEO,ALPHA,NITER,BLOCKSIZE=20) solves the reconstruction problem
   using the projection data PROJ taken over ALPHA angles, corresponding
   to the geometry descrived in GEO, using NITER iterations."""
     '\n'
     'Parameters \n'
     '-------------------------------------------------------------------\n'
     '\n'
     'proj:         Data input in the form of 3d, np.array(dtype=float32)\n'
     '\n'
     'geo:          Geometry of detector and image (see examples/Demo code)\n'
     '\n'
     'alpha:        1d array of angles corresponding to image data, np.array(dtype=float32)\n'
     '\n'
     'niter:        number of iterations of algorithm\n'
     '\n'
     'blocksize'
     'lmbda:        Sets the value of the hyperparameter.\n '
     '              Default is 1 \n'
     '\n'
     'lmbda_red:    Reduction of lambda.Every iteration \n'
     '              lambda=lambdared*lambda. Default is 0.99\n'
     '\n'
     'Init:         Describes diferent initialization techniques.\n'
     '                "none"     : Initializes the image to zeros (default) \n'
     '                "FDK"      : intializes image to FDK reconstrucition \n'
     '                "multigrid": Initializes image by solving the problem in\n'
     '                             small scale and increasing it when relative\n'
     '                             convergence is reached.\n'
     '                "image"    : Initialization using a user specified\n'
     '                             image. Not recomended unless you really\n'
     '                             know what you are doing.\n'
     '\n'
     'InitImg:      Image for the "image" initialization. (avoid)\n'
     '\n'
     'Verbose:      Feedback print statements for algorithm progress \n'
     '              default=True \n'
     '\n'
     'QualMeas:     Asks the algorithm for a set of quality measurement\n'
     '              parameters. Input should contain a list or tuple of strings of\n'
     '              quality measurement names. Example: ["CC","RMSE","MSSIM"]\n'
     '              These will be computed in each iteration.\n'
     '\n'
     'OrderStrategy:  Chooses the subset ordering strategy. Options are:\n'
     '                   "ordered"        : uses them in the input order, but divided\n'
     '                   "random"         : orders them randomply\n'
     '                   "angularDistance": chooses the next subset with the\n'
     '                                      biggest angular distance with the ones used\n'
     'Examples \n'
     '---------------------------------------------------------------------------\n'
     'See Demos/Example code for further instructions.\n'
     '---------------------------------------------------------------------------'
     '\n'
     """This file is part of the TIGRE Toolbox

        Copyright (c) 2015, University of Bath and
                            CERN-European Organization for Nuclear Research
                            All rights reserved.

        License:            Open Source under BSD.
                            See the full license at
                            https://github.com/CERN/TIGRE/license.txt

        Contact:            tigre.toolbox@gmail.com
        Codes:              https://github.com/CERN/TIGRE/
      --------------------------------------------------------------------------
        Coded by:          MATLAB (original code): Ander Biguri
                           PYTHON : Reuben Lindroos,Sam Loescher, """)
    if verbose:
        print('OS_SART algorithm in progress.')

    angleblocks, angle_index = order_subsets(alpha, blocksize, OrderStrategy)
    #     Projection weight:
    #       - fixing the geometry
    #       - making sure there are no infs in W
    geox = copy.deepcopy(geo)
    geox.sVoxel[0:] = geo.DSD - geo.DSO
    geox.sVoxel[2] = max(geox.sDetector[1], geox.sVoxel[2])
    geox.nVoxel = np.array([2, 2, 2])
    geox.dVoxel = geox.sVoxel / geox.nVoxel
    W = Ax(np.ones(geox.nVoxel, dtype=np.float32), geox, alpha, "ray-voxel")
    W[W < min(geo.dVoxel / 4)] = np.inf
    W = 1 / W

    geox = None
    #     Back_Proj weight
    #     NOTE: hstack(array,last value) as np.arange does not include upper limit of interval.
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
        xx = np.expand_dims(xx, axis=2)
        yy = np.expand_dims(yy, axis=2)

        A = (alpha + np.pi / 2)
        V = (geo.DSO / (geo.DSO + (yy * np.sin(-A)) - (xx * np.cos(-A)))) ** 2
        V = np.array(V, dtype=np.float32)


    elif geo.mode is 'parallel':
        V = np.ones([geo.nVoxel[1], geo.nVoxel[0],len(alpha)], dtype=np.float32)

    # Set up init parameters
    lq = []
    l2l =[]
    if init == 'multigrid':
        if verbose:
            print('init multigrid in progress...')
            print('default blocksize=1 for init_multigrid(OS_SART)')
        res = init_multigrid(proj, geo, alpha, alg='SART')
        if verbose:
            print('init multigrid complete.')
    if init == 'FDK':
        res=FDK(proj,geo,alpha).transpose()

    if type(init) == np.ndarray:
        if (geo.nVoxel == init.shape).all():

            res = init

        else:
            raise ValueError('wrong dimension of array for initialisation')
    elif init is None:
        res = np.zeros(geo.nVoxel, dtype=np.float32)

    # Iterate

    tic = None
    toc = time.clock()
    for i in range(niter):
        if Quameasopts is not None:
            res_prev = res
        if verbose:
            if i == 1:
                if tic is None:
                    pass
                else:
                    print('Esitmated time until completetion (s): ' + str(niter * (tic - toc)))

        for j in range(len(angleblocks)):
            if blocksize == 1:
                angle = np.array([angleblocks[j]], dtype=np.float32)
                sumax = 0
                dim_exp=True
            else:
                angle = angleblocks[j]
                sumax = 3
                dim_exp=False

            # PRESENT FOR LOOP
            res += lmbda * 1 / np.array(np.sum(np.expand_dims(
                V[:,:,angle_index[j]], axis=0), axis=sumax), dtype=np.float32) * Atb(
                adddim(W[:,:,angle_index[j]], dim_exp) * (adddim(proj[:,:,angle_index[j]], dim_exp)
                                     - Ax(res, geo, angle, 'ray-voxel')),
                geo, angle, 'FDK')
            if noneg:
                res = res.clip(min=0)

            # VERBOSE:
            # proj_err = proj[angle_index[j]] - Ax(res, geo, angle_blocks[j],'ray-voxel')
            # weighted_err = W[angle_index[j]]*proj_err
            # backprj = Atb(weighted_err,geo,angle_blocks[j],'FDK')
            # weighted_backprj = 1/V[angle_index[j]]*backprj
            # res+=1*weighted_backprj
            # res[res < 0] = 0

        if Quameasopts is not None:
            lq.append(MQ(res, res_prev, Quameasopts))
            res_prev = res
        if computel2:
            # compute l2 borm for b-Ax
            errornow=im3DNORM(proj-Ax(res,geo,alpha,'ray-voxel'),2)
            l2l.append(errornow)

        tic = time.clock()
    lmbda *= lmbda_red
    # parkerweight(projsirt,TIGRE_parameters,angles,q=1)
    if computel2:
        return res,l2l
    if Quameasopts is not None:
        return res, lq
    else:
        return res
def adddim(array,dimexp):
    # This function makes sure the dimensions of the arrays are only expanded if
    # blocksize ==1. There may be a nicer way of doing this!
    if dimexp:
        return np.expand_dims(array,axis=2)
    else:
        return array
