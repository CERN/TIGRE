#!/usr/bin/env python3
# *****************************************************************************
# *
import sys
import os
import mrcfile

import numpy as np

import scipy as sp
import tigre
from tigre.utilities import gpu, im_3d_denoise


import argparse

DESCRIPTION = 'The program will make denoising using total variation'

EXAMPLES = '''Examples:
  xmipp_tomogram_reconstruction --tiltseries ts.mrc --angles angles.xmd --thickness 300 --gpu 0 -o tomogram.mrc --method wbp --filter ram_lak
  xmipp_tomogram_reconstruction --tiltseries ts.mrc --angles angles.xmd --thickness 300 --gpu 0 -o tomogram.mrc --method fdk --filter hamming
'''


class TigreDenoisingTV:

    def __init__(self):
        # Global parameters
        self.xdim = None
        self.ydim = None
                
    @staticmethod
    def createParser():
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
            epilog=EXAMPLES,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        # --------------------- REQUIRED ARGUMENTS -----------------------------
        required = parser.add_argument_group('required arguments')
        
        required.add_argument('-i', required=True, metavar='<fnTomo>',
                              help='Volume or tomogram to be denoised')
        
        required.add_argument('-o', required=True, metavar='<fnOut>',
                              help='Output filename of the volume/tomogram.')

        # --------------------- OPTIONAL ARGUMENTS -----------------------------
        optional = parser.add_argument_group('optional arguments')
        
        optional.add_argument('--iters', type=int, metavar='<iterations>',
                              help='Number of iterations for the reconstruction algorithm.')
        
        # Nota: Se usa dest='lmbda' porque 'lambda' es una palabra reservada en Python
        optional.add_argument('--lambda', dest='lmbda', type=float, default=1.0, metavar='<lmbda>',
                              help='Hyperparameter. The update will be multiplied by this number every iteration, to make the steps bigger or smaller. Default: lmbda=1.0 for all algorithms except for irn-tv-cgls for which lmbda=5.0')
        
        optional.add_argument('--gpu', type=int, default=0, metavar='<gpuId>',
                              help='GPU Ids to be use in the image processing. (by default gpu 0) If this parameter is not set, the gpu 0 will be used')

        return parser

    def readInputParams(self):
        '''
        In this function the parameters are read. For for information about their use see the help
        '''

        self.fnTomo        = args.i
        self.fnOut         = args.o
        self.gpuId         = args.gpu if args.gpu else 0
        
        self.iterations    = args.iters if args.iters else 50
        self.lmbda         = args.lmbda if args.lmbda else 15.0

    def run(self):
        print('Starting ...')
        
        self.readInputParams()
        self.tigreDenoisingTV()

    def getGPUs(self):
        return str(self.gpuId)

    
    def run(self):
        self.readInputParams()
        
        gpuids = gpu.GpuIds()
        gpuids.devices = [int(self.getGPUs())]

        tomo = mrcfile.read(self.fnTomo)

        denoisedTomo = im_3d_denoise.im3ddenoise(np.array(tomo, dtype=np.float32), iter=self.iterations, lmbda=self.lmbda, gpuids=None)

        with mrcfile.new(self.fnOut, overwrite=True) as mrc:
            mrc.set_data(denoisedTomo)


if __name__ == '__main__':

    parser = TigreDenoisingTV.createParser()
    args   = parser.parse_args()
    TigreDenoisingTV().run()
    sys.exit(0)
