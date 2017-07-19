from __future__ import (division,print_function)
import numpy as np
from _tvdenoising import tvdenoise

def im3ddenoise(img,type,**kwargs):
    """
    IMDENOISE3D removes noise of image with different methods
       Currentyl only TV is supported. INput arguments are the iamge, the type
       of denoising ('TV' only now) and the parameters for the denoising,
       being number of iterations and hyperparameter currently available.
    --------------------------------------------------------------------------
    This file is part of the TIGRE Toolbox

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
                       PYTHON : Sam Loescher, Reuben Lindroos
    """

    tv_opts=kwargs['tv_opts']
    if type is 'TV':
        immin=min(img.ravel())
        img=(img-immin)
        immax=np.percentile(img,99)
        img=img/immax

        img0=tvdenoise(img.T,tv_opts['iter'],tv_opts['hyper'])
        img0*=immax
        img0+=immin
        return img0.T
    else:
        raise ValueError('No other options for denoising implemented')


