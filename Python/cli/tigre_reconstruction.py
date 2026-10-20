# *****************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# *****************************************************************************

import sys
import os
import argparse
import mrcfile
from utils.transforms import applyTorchTransforms
from emtable import Table
from utils.io_utils import readTltFile, readTiltSeries, readXf

import numpy as np
import torch
import tigre
from tigre.utilities.im3Dnorm import im3DNORM
import tigre.algorithms as algs
from tigre.utilities import gpu


DESCRIPTION = """\
Reconstruction of tomograms from a tilt series using multiple algorithms.

Reconstruction algorithms (--method):

  F1 - Exact algorithms:
    wbp / fbp     Weighted/Filtered Back Projection (default)
    fdk           Feldkamp-Davis-Kress

  F2 - Gradient-based:
    sirt          Simultaneous Iterative Reconstruction Technique
    sart          Simultaneous Algebraic Reconstruction Technique
    ossart        Ordered Subset SART
    asd-pocs      Adaptive Steepest Descent - POCS
    os-asd-pocs   Ordered Subset ASD-POCS
    awasd-pocs    Adaptive Weighted ASD-POCS
    pcsd          Projection-Controlled Steepest Descent
    awpcsd        Adaptive Weighted PCSD

  F3 - Krylov subspace:
    cgls          Conjugate Gradient Least Squares
    lsqr          Least Squares QR
    lsmr          Least Squares MINRES
    hybrid-lsqr   Hybrid LSQR
    ab-gmres      AB-Generalized Minimal Residual
    ba-gmres      BA-Generalized Minimal Residual
    irn-tv-cgls   IRN TV CGLS

  F4 - Statistical:
    mlem          Maximum Likelihood Expectation Maximization

  F5 - Variational:
    fista         Fast Iterative Shrinkage Thresholding
    sart-tv       SART with Total Variation
"""

EXAMPLES = """\
Examples:
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method wbp --filter ram_lak
  %(prog)s --tiltseries ts.mrc --angles angles.xmd --thickness 300 --gpu 0 -o tomogram.mrc --method fdk --filter hamming
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method sirt --iter 20
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method sart
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method ossart
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method lsmr
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method asd-pocs
  %(prog)s --tiltseries ts.mrc --angles angles.tlt --thickness 300 --gpu 0 -o tomogram.mrc --method cgls
"""

class TomogramReconstruction:

    def __init__(self):
        self.fnTs          = None
        self.fnAngles      = None
        self.thickness     = None
        self.method        = None
        self.fnOut         = None
        self.normalizeTi   = ''
        self.backprojector = None
        self.useTigreInterpolation = False
        self.gpuId         = 0

        self.xdim = None
        self.ydim = None


    @staticmethod
    def createParser():
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
            epilog=EXAMPLES,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        required = parser.add_argument_group('required arguments')
        required.add_argument('--tiltseries', required=True, metavar='FILE',
                              help='Tilt series file (.mrc, .mrcs, .ali, or .st)')
        required.add_argument('--angles', required=True, metavar='FILE',
                              help='Angles file (.tlt) with tilt angles')
        required.add_argument('--xf', required=True, metavar='FILE',
                              help='text file (.xf) with in plane rotations (as rotation matrix) and shifts')
        required.add_argument('--thickness', required=True, type=int, metavar='N',
                              help='Thickness in pixels of the reconstructed tomogram')
        required.add_argument('-o', required=True, metavar='FILE',
                              help='Output filename of the tomogram')

        parser.add_argument('--method', default='wbp', metavar='METHOD',
                            help='Reconstruction algorithm (default: wbp). See description for full list.')
        parser.add_argument('--filter', default=None, metavar='FILTER',
                            help='Filter for FDK or FBP: ram_lak (default), shepp_logan, cosine, hamming, hann')
        parser.add_argument('--notNormalize', action='store_true', 
                            help='Normalisation mode. Use "standard" for zero mean and unit std dev.')
        parser.add_argument('--tigreInterpolation', action='store_true',
                            help='Set this flag if Tigre should  manage the tilt series alignment')
        parser.add_argument('--gpu', type=int, default=0, metavar='ID',
                            help='GPU id to use (default: 0)')

        iterative = parser.add_argument_group('iterative algorithm options')
        iterative.add_argument('--iter', type=int, default=None, metavar='N',
                               help='Number of iterations (default depends on method)')
        iterative.add_argument('--lambda', dest='lmbda', type=float, default=None, metavar='F',
                               help='Step size hyperparameter (default: 1.0, or 5.0 for irn-tv-cgls)')
        iterative.add_argument('--lambdared', type=float, default=None, metavar='F',
                               help='Lambda reduction multiplier per iteration (default: 0.999)')
        iterative.add_argument('--blocksize', type=int, default=None, metavar='N',
                               help='Projection block size for ossart/os-asd-pocs/awasd-pocs (default: 10-20)')
        iterative.add_argument('--order', default=None, metavar='STRATEGY',
                               help='Subset ordering for ossart: ordered, random, angularDistance (default: random)')
        iterative.add_argument('--backprojector', default=None, metavar='NAME',
                               help='Backprojector for ab-gmres / ba-gmres')

        tv = parser.add_argument_group('TV-based algorithm options')
        tv.add_argument('--alpha', type=float, default=None, metavar='F',
                        help='TV hyperparameter (default: 0.002)')
        tv.add_argument('--tviter', type=int, default=None, metavar='N',
                        help='TV denoising iterations (default depends on method)')
        tv.add_argument('--tvlambda', type=float, default=None, metavar='F',
                        help='TV multiplier (default depends on method)')
        tv.add_argument('--alpha_red', type=float, default=None, metavar='F',
                        help='TV hyperparameter reduction rate (default: 0.95)')
        tv.add_argument('--ratio', type=float, default=None, metavar='F',
                        help='Max allowed image/TV update ratio (default: 0.94)')

        return parser


    def readInputParams(self, args):
        self.fnTs        = args.tiltseries
        self.fnAngles    = args.angles
        self.fnXf        = args.xf
        self.thickness   = args.thickness
        self.method      = self.fixMethod(args.method.lower())
        self.fnOut       = args.o
        self.gpuId       = args.gpu
        self.normalizeTi = not args.notNormalize
        self.useTigreInterpolation = args.tigreInterpolation

        if self.method == 'wbp':
            self.method = 'fbp'

        if self.method in ('fdk', 'fbp'):
            filterArg = args.filter if args.filter else 'ram_lak'
            self.checkFilter(filterArg.lower(), self.method)

        if self.method in ('sirt', 'sart', 'ossart'):
            self.iterations = args.iter      if args.iter      is not None else 20
            self.lmbda      = args.lmbda     if args.lmbda     is not None else 1.0
            self.lambdared  = args.lambdared if args.lambdared is not None else 0.999
            self.qualmeas   = ["RMSE", "SSD"]

        if self.method == 'ossart':
            self.blocksize = args.blocksize if args.blocksize is not None else 10
            self.order     = args.order     if args.order     is not None else "random"

        if self.method in ('cgls', 'lsqr', 'lsmr', 'hybridlsqr', 'abgmres', 'bagmres'):
            self.iterations = args.iter if args.iter is not None else 30

        if self.method in ('abgmres', 'bagmres'):
            self.backprojector = args.backprojector

        if self.method == 'asdpocs':
            self.iterations = args.iter if args.iter is not None else 30

        if self.method in ('lsmr', 'asdpocs', 'osasdpocs', 'awasdpocs'):
            self.lmbda = args.lmbda if args.lmbda is not None else 1.0

        if self.method in ('asdpocs', 'osasdpocs', 'awasdpocs'):
            self.alpha     = args.alpha     if args.alpha     is not None else 0.002
            self.tviter    = args.tviter    if args.tviter    is not None else 25
            self.lambdared = args.lambdared if args.lambdared is not None else 0.9999
            self.alpha_red = args.alpha_red if args.alpha_red is not None else 0.95
            self.ratio     = args.ratio     if args.ratio     is not None else 0.94

        if self.method in ('pcsd', 'awpcsd'):
            self.iterations = args.iter if args.iter is not None else 20

        if self.method == 'osasdpocs':
            self.blocksize  = args.blocksize if args.blocksize is not None else 10
            self.iterations = args.iter      if args.iter      is not None else 20

        if self.method == 'awasdpocs':
            self.blocksize = args.blocksize if args.blocksize is not None else 20
            self.iterations = args.iter      if args.iter      is not None else 20


        if self.method == 'irntvcgls':
            self.lmbda       = args.lmbda if args.lmbda is not None else 5.0
            self.iterations  = args.iter  if args.iter  is not None else 10
            self.niter_outer = 2

        if self.method == 'fista':
            self.iterations = args.iter     if args.iter     is not None else 100
            self.tviter     = args.tviter   if args.tviter   is not None else 100
            self.tvlambda   = args.tvlambda if args.tvlambda is not None else 20

        if self.method == 'sarttv':
            self.iterations = args.iter     if args.iter     is not None else 30
            self.tviter     = args.tviter   if args.tviter   is not None else 50
            self.tvlambda   = args.tvlambda if args.tvlambda is not None else 50
            self.alpha_red  = args.alpha_red if args.alpha_red is not None else 0.95

        if self.method == 'mlem':
            self.iterations = args.iter if args.iter is not None else 500
        


    def fixMethod(self, method):
        method = method.replace("-", "")
        method = method.replace("_", "")
        return method


    def fixFilter(self, candidate, method):
        listFilters = ["ram_lak", "shepp_logan", "cosine", "hamming", "hann"]
        if '-' in candidate:
            candidate = candidate.replace('-', '_')
        if not (candidate in listFilters):
            raise Exception('The selected filter presents a problem, please check the --filter')
        return candidate, listFilters


    def checkFilter(self, candidate, method):
        candidate, listFilters = self.fixFilter(candidate, method)
        if candidate in listFilters:
            self.filterToApply = listFilters[listFilters.index(candidate)]
        else:
            raise Exception('The selected filter does not exist, please check the --filter flag')
        print(self.filterToApply)



    def getGPUs(self):
        return str(self.gpuId)


    def run(self, args):
        print('Starting ...')
        from scipy.spatial.transform import Rotation as R

        self.readInputParams(args)
        ts, self.xdim, self.ydim = readTiltSeries(self.fnTs, normalizeTis=self.normalizeTi)

        device = self.gpuId or ("cuda" if torch.cuda.is_available() else "cpu")

        tiltAngles = []
        rotAngles  = []
        offSets    = []

        if os.path.splitext(self.fnAngles)[1] == '.tlt':
            tiltAngles = readTltFile(self.fnAngles)

            transforms = readXf(self.fnXf)
            if len(transforms) != ts.shape[0]:
                raise ValueError(
                    f"Number of transformations ({len(transforms)}) does not match "
                    f"number of images ({ts.shape[0]})"
                )

            matrices2x2 = np.stack([m for m, _ in transforms])
            rotAngles = np.stack([np.arctan2(m2d[0, 0], -m2d[0, 1]) for m2d in matrices2x2])
            offSets = np.stack([t for _, t in transforms])

        
        tiltAngles = np.array(tiltAngles) * np.pi / 180.0
        
        if self.useTigreInterpolation:
            print('Tigre will manage the tilt series alignment')
            offSets = -offSets
            print(offSets)
            rotAngles = -rotAngles
            pass
        else:
            print('Applying the tilt series alignment')
            ts = applyTorchTransforms(ts, matrices2x2, offSets, device=device)
            rotAngles = None
            offSets = None
        
        ts = ts.astype(np.float32)
        # with mrcfile.new('interT.mrcs', overwrite=True) as mrc:
        #     mrc.set_data(ts2)
        
        self.tigreReconstruction(ts, tiltAngles, rotAngles=rotAngles, offSets=offSets, recMethod=self.method)


    def tigreReconstruction(self, ts, tiltAngles, rotAngles=None, offSets=None, recMethod='fbp'):


        gpuids = gpu.GpuIds()
        gpuids.devices = [int(self.getGPUs())]

        geo = tigre.geometry(mode="parallel", nVoxel=np.array([self.ydim, self.xdim, self.thickness]))
        '''
        if self.useTigreInterpolation:
            if rotAngles is not None:
                geo.rotDetector = rotAngles
            if offSets is not None:
                geo.offDetector = offSets
        '''

        if recMethod == 'fbp':
            reconstruction = algs.fbp(ts, geo, tiltAngles, filter=self.filterToApply, noneg=False, gpuids=gpuids)

        elif recMethod == 'fdk':
            reconstruction = algs.fdk(ts, geo, tiltAngles, filter=self.filterToApply, noneg=False, gpuids=gpuids)

        elif recMethod == 'sirt':
            reconstruction = algs.sirt(ts, geo, tiltAngles, self.iterations, lmbda=self.lmbda, lmbda_red=self.lambdared, verbose=True, noneg=False, gpuids=gpuids)

        elif recMethod == 'sart':
            reconstruction = algs.sart(ts, geo, tiltAngles, self.iterations, lmbda=self.lmbda, lmbda_red=self.lambdared, verbose=True, noneg=False, gpuids=gpuids)

        elif recMethod == 'ossart':
            reconstruction = algs.ossart(ts, geo, tiltAngles, self.iterations, lmbda=self.lmbda, lmbda_red=self.lambdared, verbose=True, noneg=False, blocksize=self.blocksize, OrderStrategy=self.order, gpuids=gpuids)

        elif recMethod == 'cgls':
            reconstruction = algs.cgls(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)

        elif recMethod == 'lsqr':
            reconstruction = algs.lsqr(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)

        elif recMethod == 'lsmr':
            reconstruction = algs.lsmr(ts, geo, tiltAngles, self.iterations, lmbda=self.lmbda, noneg=False, gpuids=gpuids)

        elif recMethod == 'hybridlsqr':
            reconstruction = algs.hybrid_lsqr(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)

        elif recMethod == 'abgmres':
            if self.backprojector is None:
                reconstruction = algs.ab_gmres(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)
            else:
                reconstruction = algs.ab_gmres(ts, geo, tiltAngles, self.iterations, backprojector="FDK", noneg=False, gpuids=gpuids)

        elif recMethod == 'bagmres':
            if self.backprojector is None:
                reconstruction = algs.ba_gmres(ts, geo, tiltAngles, self.iterations, noneg=False)
            else:
                reconstruction = algs.ba_gmres(ts, geo, tiltAngles, self.iterations, backprojector="FDK", noneg=False, gpuids=gpuids)

        elif recMethod == 'asdpocs':
            epsilon = im3DNORM(tigre.Ax(algs.fdk(ts, geo, tiltAngles), geo, tiltAngles) - ts, 2) * 0.15
            reconstruction = algs.asd_pocs(ts, geo, tiltAngles, self.iterations, tviter=self.tviter, maxl2err=epsilon,
                                            alpha=self.alpha, lmbda=self.lmbda, lmbda_red=self.lambdared, rmax=self.ratio, verbose=True, noneg=False, gpuids=gpuids)

        elif recMethod == 'osasdpocs':
            epsilon = im3DNORM(tigre.Ax(algs.fdk(ts, geo, tiltAngles), geo, tiltAngles) - ts, 2) * 0.15
            reconstruction = algs.os_asd_pocs(ts, geo, tiltAngles, self.iterations, tviter=self.tviter, maxl2err=epsilon,
                                               alpha=self.alpha, lmbda=self.lmbda, lmbda_red=self.lambdared, rmax=self.ratio,
                                               verbose=True, blocksize=self.blocksize, noneg=False, gpuids=gpuids)

        elif recMethod == 'awasdpocs':
            epsilon = im3DNORM(tigre.Ax(algs.fdk(ts, geo, tiltAngles), geo, tiltAngles) - ts, 2) * 0.15
            reconstruction = algs.awasd_pocs(ts, geo, tiltAngles, self.iterations, tviter=self.tviter, maxl2err=epsilon,
                                              alpha=self.alpha, lmbda=self.lmbda, lmbda_red=self.lambdared, rmax=self.ratio,
                                              verbose=True, delta=np.array([-0.005]), noneg=False, gpuids=gpuids)

        elif recMethod == 'irntvcgls':
            reconstruction = algs.irn_tv_cgls(ts, geo, tiltAngles, self.iterations, lmbda=self.lmbda, niter_outer=self.niter_outer, noneg=False, gpuids=gpuids)

        elif recMethod == 'fista':
            reconstruction = algs.fista(ts, geo, tiltAngles, self.iterations, tviter=self.tviter, tvlambda=self.tvlambda, noneg=False, gpuids=gpuids)

        elif recMethod == 'sarttv':
            reconstruction = algs.sart_tv(ts, geo, tiltAngles, self.iterations, tvlambda=self.tvlambda, tviter=self.tviter, noneg=False, gpuids=gpuids)

        elif recMethod == 'mlem':
            reconstruction = algs.mlem(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)

        elif recMethod == 'pcsd':
            reconstruction = algs.pcsd(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)

        elif recMethod == 'awpcsd':
            reconstruction = algs.aw_pcsd(ts, geo, tiltAngles, self.iterations, noneg=False, gpuids=gpuids)

        else:
            raise Exception('The selected reconstruction method does not exist, check the --method flag')

        with mrcfile.new(self.fnOut, overwrite=True) as mrc:
            mrc.set_data(np.transpose(reconstruction, (2, 0, 1)))


if __name__ == '__main__':
    parser = TomogramReconstruction.createParser()
    args   = parser.parse_args()
    rec    = TomogramReconstruction()
    rec.run(args)
    sys.exit(0)
