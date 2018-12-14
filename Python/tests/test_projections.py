import unittest
from tigre.utilities.geometry import geometry
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
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
        self.img = np.ones((63, 62, 61), dtype=np.float32)
        self.proj = Ax(self.img,self.geo,self.angles)
    def test_geo(self):
        self.assertIsNone(self.geo.check_geo(self.angles))

    def test_proj_interpolated_shape(self):
        img = self.img
        self.assertTupleEqual(Ax(img,self.geo,self.angles,'interpolated').shape,(100,128,127))

    def test_proj_ray_voxel_shape(self):
        self.assertTupleEqual(Ax(self.img,self.geo,self.angles,'ray-voxel').shape,(100,128,127))

    def test_proj_parallel_interpolated_shape(self):
        setattr(self.geo,'mode','parallel')
        self.assertTupleEqual(Ax(self.img,self.geo,self.angles,'interpolated').shape,(100,128,127))

    def test_proj_parallel_ray_voxel_shape(self):
        setattr(self.geo,'mode','parallel')
        self.assertTupleEqual(Ax(self.img,self.geo,self.angles,'ray-voxel').shape,(100,128,127))

    def test_backproj_matched_cone_shape(self):
        setattr(self.geo,'mode','cone')
        self.assertTupleEqual(Atb(self.proj,self.geo,self.angles,'matched').shape,(63,62,61))

    def test_backproj_FDK_cone_shape(self):
        setattr(self.geo,'mode','cone')
        self.assertTupleEqual(Atb(self.proj, self.geo, self.angles, 'FDK').shape, (63, 62, 61))

    def test_backproj_parallel_shape(self):
        setattr(self.geo,'mode','parallel')
        self.assertTupleEqual(Atb(self.proj, self.geo, self.angles,'matched').shape, (63, 62, 61))

    def test_backproj_matched_cone_spherical_shape(self):
        setattr(self.geo,'mode','cone')
        angles_1 = self.angles
        angles_2 = np.ones((len(angles_1)),dtype=np.float32)*np.array([np.pi/4],dtype=np.float32)
        angles_3 = np.zeros((len(angles_1)),dtype=np.float32)
        new_angles = np.vstack((angles_1,angles_2,angles_3)).T
        self.assertTupleEqual(Atb(self.proj, self.geo, new_angles,'matched').shape, (63, 62, 61))

    def test_backproj_FDK_cone_spherical_shape(self):
        setattr(self.geo,'mode','cone')
        angles_1 = self.angles
        angles_2 = np.ones((len(angles_1)), dtype=np.float32) * np.array([np.pi / 4], dtype=np.float32)
        angles_3 = np.zeros((len(angles_1)), dtype=np.float32)
        new_angles = np.vstack((angles_1, angles_2, angles_3)).T
        self.assertTupleEqual(Atb(self.proj, self.geo, new_angles,'FDK').shape, (63, 62, 61))

    def test_backproj_parallel_spherical_shape(self):
        setattr(self.geo,'mode','parallel')
        angles_1 = self.angles
        angles_2 = np.ones((len(angles_1)), dtype=np.float32) * np.array([np.pi / 4], dtype=np.float32)
        angles_3 = np.zeros((len(angles_1)), dtype=np.float32)
        new_angles = np.vstack((angles_1, angles_2, angles_3)).T
        self.assertTupleEqual(Atb(self.proj, self.geo, new_angles).shape, (63, 62, 61))

if __name__=='__main__':
        unittest.main()
