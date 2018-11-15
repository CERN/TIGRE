import unittest
from tigre.geometry import geometry
from tigre.Ax import Ax
from tigre.Atb import Atb
import numpy as np

class AlgorithmsTestCase(unittest.TestCase):

    def setUp(self):
        geo = geometry()
        geo.DSD = 1536
        geo.DSO = 1000
        geo.nDetector = np.array((128, 127))
        geo.dDetector = np.array((0.8, 0.8)) * 4
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
    def test_geo(self):
        self.assertIsNone(self.geo.check_geo(self.angles))
