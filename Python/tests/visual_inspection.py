from __future__ import print_function
import os
import sys

import tigre
import numpy as np
import tigre.demos.Test_data.data_loader as data_loader
import tigre.algorithms as algs
from tigre.utilities.Ax import Ax
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
from matplotlib import pyplot as plt

# ---------------------PLOT---------------------------
def plot_algs(alglist, proj, geo, angles, niter=10, **kwargs):
    l= []
    for alg in alglist:
        if alg in ['fbp','fdk','FDK']:
            l.append(getattr(algs, alg)(proj, geo, angles)[geo.nVoxel[0]/2])
        else:
            kwargs.update(dict(sup_kw_warning=True))
            l.append(getattr(algs, alg)(proj, geo, angles, niter, **kwargs)[geo.nVoxel[0]/2])
    for i in range(len(l)):
        plt.matshow(l[i])
        plt.title(str(alglist[i]) + ' ' + str(geo.mode))
        if kwargs.has_key('colorbar'):
            plt.colorbar()
    plt.show()


def input_parser(sysarglist):

    if len(sysarglist) < 2:
        raise ValueError('visual test requires at least 2 input parameters, type --help')
    if sys.argv[1]== '--help':
        print('arg 1: mode: cone, parallel \n' +
              'arg 2: alglist, eg: "ossart sirt fdk" \n'
              'arg 3: kwargs, eg: "blocksize=20, nVoxel=32,64,128')
        raise SystemExit()

    kwargs = dict(niter=10,
                  nVoxel=np.array([64, 64, 64]),
                  colorbar=True,
                  krylov='interpolated')

    kwargs.update(dict(mode=sysarglist[1],
                       alglist=sysarglist[2]))
    if len(sysarglist) > 3:
        for arg in sysarglist[3].split(' '):
            k = arg.split("=")[0]
            v = arg.split("=")[1]
            if k == 'nVoxel':
                v = np.array([int(val) for val in v.split(',')])
            kwargs.update({k: v})
            print(kwargs)

    return kwargs


def setUp(**kwargs):
    if kwargs.get('mode') == 'cone':
        geo = tigre.geometry(mode='cone', nVoxel=kwargs.get('nVoxel'), default=True)
        angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
        angles_3 = np.zeros((100), dtype=np.float32)
        angles = np.vstack((angles_1, angles_3, angles_3)).T
        img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
        proj = tigre.Ax(img, geo, angles, kwargs.get('krylov'))
        if kwargs.get('alglist') == 'all':
            alglist = ['sart', 'sirt', 'ossart', 'asd_pocs', 'FDK', 'cgls']
        else:
            alglist = kwargs.pop('alglist').split()

        plot_algs(alglist, proj, geo, angles, niter=int(kwargs.pop('niter')), **kwargs)

    if kwargs.get('mode') == 'parallel':
        geo = tigre.geometry(mode='parallel', nVoxel=kwargs.get('nVoxel'), default=True)
        angles_1 = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
        angles_3 = np.zeros((100), dtype=np.float32)
        angles = np.vstack((angles_1, angles_3, angles_3)).T
        img = data_loader.load_head_phantom(number_of_voxels=geo.nVoxel)
        proj = tigre.Ax(img, geo, angles, kwargs.get('krylov'))
        if kwargs.get('alglist') == 'all':
            alglist = ['sart', 'sirt', 'ossart', 'asd_pocs', 'fbp', 'cgls']
        else:
            alglist = kwargs.pop('alglist').split()

        plot_algs(alglist, proj, geo, angles, niter=int(kwargs.pop('niter')), **kwargs)


if __name__ == '__main__':
    setUp(**input_parser(sys.argv))
