from __future__ import division
import unittest
import tigre
from tigre.utilities.Ax import Ax
import numpy as np
from tigre.algorithms.iterative_recon_alg import IterativeReconAlg
import tigre.algorithms as algs
class AlgorithmsTestCase(unittest.TestCase):

    def setUp(self):
        geo = tigre.geometry()
        geo.DSD = 1536
        geo.DSO = 1000
        geo.nDetector = np.array((128, 127))
        geo.dDetector = np.array((0.8, 0.8)) * 4.
        geo.sDetector = geo.nDetector * geo.dDetector
        geo.nVoxel = np.array((63, 62, 61))
        geo.sVoxel = np.array((256, 256, 256))
        geo.dVoxel = geo.sVoxel / geo.nVoxel
        geo.offOrigin = np.array((0, 0, 0))
        geo.offDetector = np.array((0, 0))
        geo.accuracy = 0.5
        geo.mode = 'cone'
        self.angles = np.linspace(0, 2 * np.pi, 100, dtype=np.float32)
        self.geo = geo
        self.img = np.ones((63, 62, 61), dtype=np.float32)
        self.proj = Ax(self.img,self.geo,self.angles)
        self.niter=2
    def test_mainiter_itrecalg(self):
        alg = IterativeReconAlg(self.proj,self.geo,self.angles,self.niter,**dict(verbose=False))
        self.assertIsNone(alg.run_main_iter())

    def test_shape_sart(self):
        self.assertTupleEqual(tigre.algorithms.sart(self.proj,self.geo,self.angles,niter=1).shape,
                              tuple(self.geo.nVoxel))
    def test_shape_sirt(self):
        self.assertTupleEqual(tigre.algorithms.sirt(self.proj,self.geo,self.angles,niter=1).shape,
                              tuple(self.geo.nVoxel))
    def test_shape_ossart(self):
        self.assertTupleEqual(tigre.algorithms.sirt(self.proj,self.geo,self.angles,niter=1).shape,
                              tuple(self.geo.nVoxel))
    def test_shape_asdpocs(self):
        self.assertTupleEqual(algs.asd_pocs(self.proj,self.geo,self.angles,niter=1).shape,
                              tuple(self.geo.nVoxel))
    def test_shape_awasdpocs(self):
        self.assertTupleEqual(algs.awasd_pocs(self.proj,self.geo,self.angles,niter=1).shape,
                              tuple(self.geo.nVoxel))
if __name__=='__main__':
        unittest.main()