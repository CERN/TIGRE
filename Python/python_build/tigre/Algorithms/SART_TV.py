from __future__ import division
from __future__ import print_function
import time
from tigre.Utilities.init_multigrid import init_multigrid
from scipy.linalg import *
import numpy as np
import copy
from _Ax import Ax
from _Atb import Atb
from tigre.Utilities.im3Dnorm import im3DNORM
from tigre.Utilities.order_subsets import order_subsets
from tigre.Utilities.Measure_Quality import Measure_Quality as MQ
from tigre.Utilities.im3ddenoise import im3ddenoise
def SART_TV(proj, geo, alpha, niter,
         lmbda=1, lmbda_red=0.99, OrderStrategy=None, Quameasopts=None, init=None, verbose=True,noneg=True,computel2=False,**kwargs):

    ('\n'
     """SART_TV solves Cone Beam CT image reconstruction using Oriented Subsets
              Simultaneous Algebraic Reconxtruction Techique algorithm

   SART(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
   using the projection data PROJ taken over ALPHA angles, corresponding
   to the geometry described in GEO, using NITER iterations."""
     '\n'
     'Parameters \n'
     '-------------------------------------------------------------------\n'
     '\n'
     'proj:         Data input in the form of 3d np.array(dtype=float32)\n'
     '\n'
     'geo:          Geometry of detector and image (see examples/Demo code)\n'
     '\n'
     'alpha:        1d array of angles corresponding to image data\n'
     '\n'
     'niter:        number of iterations of algorithm\n'
     '\n'
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
                           PYTHON : Sam Loescher, Reuben Lindroos""")
    if verbose:
        print('SART_TV algorithm in progress.')
    d_opts={'tvopts':{'hyper':15,'iter':50}}
    d_opts.update(kwargs)
    blocksize = 1
    angles, angle_index = order_subsets(alpha, blocksize, OrderStrategy)
    alpha = angles.ravel()

    #     Projection weight:
    #       - fixing the geometry
    #       - making sure there are no infs in W
    geox = copy.deepcopy(geo)
    geox.sVoxel[0:]=geo.DSD-geo.DSO
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

        xv =np.arange(start, stop+step, step)

        start = geo.sVoxel[2] / 2 - geo.dVoxel[2] / 2 + geo.offOrigin[2]
        stop = -geo.sVoxel[2] / 2 + geo.dVoxel[2] / 2 + geo.offOrigin[2]
        step = -geo.dVoxel[2]

        yv = -1 * np.arange(start, stop+step, step)

        (yy,xx) = np.meshgrid(yv,xv)
        xx = np.expand_dims(xx,axis=2)
        yy = np.expand_dims(yy,axis=2)
        A = (alpha + np.pi / 2)
        V = (geo.DSO / (geo.DSO + (yy * np.sin(-A)) - (xx * np.cos(-A)))) ** 2
        V = np.array(V, dtype=np.float32)




    elif geo.mode == 'parallel':
        V = np.ones([geo.nVoxel[1], geo.nVoxel[0],len(alpha)], dtype=np.float32)

    # Set up init parameters
    if init == 'multigrid':
        if verbose:
            print('init multigrid in progress...')
        res = init_multigrid(proj, geo, alpha, alg='SART')
        if verbose:
            print('init multigrid complete.')
    if init == 'FDK':
        raise ValueError('FDK not implemented as of yet (coming soon)!')

    if type(init) == np.ndarray:
        if (geo.nVoxel == init.shape).all():

            res = init

        else:
            raise ValueError('wrong dimension of array for initialisation')
    elif init is None:
        res = np.zeros(geo.nVoxel, dtype=np.float32)

    # List for storing Quality measure
    lq = []
    l2l=[]
    # Iterate


    for i in range(niter):
        if Quameasopts is not None:
            res_prev = res
        if verbose:
            if i == 1:
                if tic is None:
                    pass
                else:
                    print('Esitmated time until completetion (s): ' + str((niter - 1) * (tic - toc)))

        toc = time.clock()
        for j in range(len(alpha)):

            # VERBOSE:
            # angle = np.array([alpha[j]], dtype=np.float32)
            # proj_err = proj[angle_index[j]] - Ax(res, geo, angle, 'ray-voxel')
            # weighted_err = W[angle_index[j]] * proj_err
            # backprj = Atb(weighted_err, geo, angle, 'FDK')
            # weighted_backprj = 1 / V[angle_index[j]] * backprj
            # res += weighted_backprj
            # res[res<0]=0

            angle = np.array([alpha[j]], dtype=np.float32)

            res += lmbda * 1 / V[:,:,angle_index[j]] * Atb(
                np.expand_dims(W[:,:,j],axis=2) * (np.expand_dims(proj[:,:,angle_index[j]],axis=2) - Ax(res, geo, angle, 'ray-voxel')),
                geo, angle, 'FDK')
            if noneg:
                res = res.clip(min=0)

        if Quameasopts is not None:
            lq.append(MQ(res, res_prev, Quameasopts))
        if computel2:
            # compute l2 borm for b-Ax
            errornow=im3DNORM(proj-Ax(res,geo,alpha,'ray-voxel'),2)
            l2l.append(errornow)

            res_prev = res
        tic = time.clock()
    lmbda *= lmbda_red
    res=im3ddenoise(res,'TV',tv_opts=d_opts['tvopts'])
    if computel2:
        return res,l2l
    if Quameasopts is not None:
        return res, lq
    else:
        return res