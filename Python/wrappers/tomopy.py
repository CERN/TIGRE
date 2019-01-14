from __future__ import (absolute_import, division, print_function)
from tigre.geometry import TIGREParameters
import numpy as np
import tomopy
import copy
from matplotlib import pyplot as plt
import tigre
from _Ax import Ax

# TODO: modify 'center' so that it actually does something.

def tigre_w(tomo,  theta, recon, center=None, **kwargs):
    '''
    Parameters:
    -----------
    tomo: object being reconstructed (np.array)
    theta: angles (np.array)
    recon: object being reconstructed to (np.array)
    center: default None
    kwargs: options for tigre see help(tigre.ex_func) for further instructions

    Examples:
    --------
    >>> import tomopy
    >>> import tigre
    >>> obj = tomopy.shepp3d()
    >>> ang = np.linspace(0,2*pi,100)
    >>> sim = tomopy.project(obj,ang)
    >>> rec = tigre_w(sim,ang,obj,options={'method':'SART_TV','num_iter':10,
    >>>                                    'geo':{'nDetector:'[sim.shape[1],sim.shape[0]]}})
    >>> tigre.plotImg(rec)
    >>> # returns plot of reconstructed image.
    -----------------------------------------------------------------
    Coded by:          Reuben Lindroos
    '''
    default_opts = {'blocksize': 20,
                    'lmbda': 1,
                    'lmbda_red': 0.99,
                    'OrderStrategy': None,
                    'Quameasopts': None,
                    'init': None,
                    'verbose': True,
                    'noneg': True,
                    'computel2': False,
                    'geo': {}
                    }

    opts = default_opts
    opts.update(kwargs['options'])

    niter = opts['num_iter']
    m_opt = opts['method']

    # generate tigre geometry

    geo=tigre.geo(recon.shape,geo=opts['geo'])
    
    if m_opt is 'OS_SART':
        res = tigre.OS_SART(tomo, geo, theta, niter,
                            blocksize=opts['blocksize'],
                            lmbda=opts['lmbda'],
                            lmbda_red=opts['lmbda_red'],
                            OrderStrategy=opts['OrderStrategy'],
                            Quameasopts=opts['Quameasopts'],
                            init=opts['init'],
                            verbose=opts['verbose'],
                            computel2=opts['computel2'])

        return res
    if m_opt is 'SART':
        res = tigre.SART(tomo, geo, theta, niter,
                         lmbda=opts['lmbda'],
                         lmbda_red=opts['lmbda_red'],
                         OrderStrategy=opts['OrderStrategy'],
                         Quameasopts=opts['Quameasopts'],
                         init=opts['init'],
                         verbose=opts['verbose'],
                         computel2=opts['computel2'])
        return res
    if m_opt is 'SIRT':
        res = tigre.SIRT(tomo, geo, theta, niter,
                         lmbda=opts['lmbda'],
                         lmbda_red=opts['lmbda_red'],
                         OrderStrategy=opts['OrderStrategy'],
                         Quameasopts=opts['Quameasopts'],
                         init=opts['init'],
                         verbose=opts['verbose'],
                         computel2=opts['computel2'])
    if m_opt is 'SART_TV':
        res = tigre.SART_TV(tomo, geo, theta, niter,
                             lmbda=opts['lmbda'],
                             lmbda_red=opts['lmbda_red'],
                             OrderStrategy=opts['OrderStrategy'],
                             Quameasopts=opts['Quameasopts'],
                             init=opts['init'],
                             verbose=opts['verbose'],
                             computel2=opts['computel2'])
        return res
    if m_opt is 'FDK':
        res = tigre.FDK(tomo, geo, theta)
        return res
    if m_opt is 'FBP':
        res = tigre.FBP(tomo, geo, theta)
        return res
    if m_opt not in ['SIRT', 'SART', 'OS_SART', 'FDK','SART_TV', 'FBP']:
        raise ValueError('Algorithm for TIGRE not recognised')


