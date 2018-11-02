import unittest
from tigre.geometry import geometry
from tigre.Ax import Ax
from tigre.Atb import Atb
import numpy as np

class ProjectionTestCase(unittest.TestCase):

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

    def test_geo(self):
        self.assertIsNone(self.geo.check_geo(self.angles))

    def test_proj_shape(self):
        img = np.ones((63,62,61),dtype=np.float32)
        self.assertTupleEqual(Ax(img,self.geo,self.angles,'interpolated').shape,(100,128,127))

    def test_backproj_shape(self):
        img = np.ones((63,62,61),dtype=np.float32)
        proj = Ax(img,self.geo,self.angles,'interpolated')
        self.assertTupleEqual(Atb(proj,self.geo,self.angles).shape,(63,62,61))

if __name__=='__main__':
        unittest.main()
