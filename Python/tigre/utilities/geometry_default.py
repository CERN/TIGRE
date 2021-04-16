from __future__ import division

import numpy as np
from tigre.utilities.geometry import Geometry


class ConeGeometryDefault(Geometry):
    def __init__(self, high_resolution=True, nVoxel=None):

        Geometry.__init__(self)
        if high_resolution:
            # VARIABLE                                          DESCRIPTION                    UNITS
            # -------------------------------------------------------------------------------------
            self.DSD = 1536.0  # Distance Source Detector      (mm)
            self.DSO = 1000.0  # Distance Source Origin        (mm)
            # Detector parameters
            self.nDetector = np.array((512, 512))  # number of pixels              (px)
            self.dDetector = np.array((0.8, 0.8))  # size of each pixel            (mm)
            self.sDetector = self.nDetector * self.dDetector  # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((512, 512, 512))  # number of voxels              (vx)
            self.sVoxel = np.array((256, 256, 256))  # total size of the image       (mm)
            self.dVoxel = self.sVoxel / self.nVoxel  # size of each voxel            (mm)
            # Offsets
            self.offOrigin = np.array((0, 0, 0))  # Offset of image from origin   (mm)
            self.offDetector = np.array((0, 0))  # Offset of Detector            (mm)

            # Auxiliary
            self.accuracy = 0.5  # Accuracy of FWD proj          (vx/sample)  # noqa: E501
            # Mode
            self.mode = "cone"  # parallel, cone                ...
            self.filter = None
        else:
            # VARIABLE                                          DESCRIPTION                    UNITS
            # -------------------------------------------------------------------------------------
            self.DSD = 1536.0  # Distance Source Detector      (mm)
            self.DSO = 1000.0  # Distance Source Origin        (mm)
            # Detector parameters
            self.nDetector = np.array((128, 128))  # (V,U) number of pixels        (px)
            self.dDetector = np.array((0.8, 0.8)) * 4  # size of each pixel            (mm)
            self.sDetector = self.nDetector * self.dDetector  # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((64, 64, 64))  # number of voxels              (vx)
            self.sVoxel = np.array((256, 256, 256))  # total size of the image       (mm)
            self.dVoxel = self.sVoxel / self.nVoxel  # size of each voxel            (mm)
            # Offsets
            self.offOrigin = np.array((0, 0, 0))  # Offset of image from origin   (mm)
            self.offDetector = np.array((0, 0))  # Offset of Detector            (mm)
            self.rotDetector = np.array((0, 0, 0))
            # Auxiliary
            self.accuracy = 0.5  # Accuracy of FWD proj          (vx/sample)  # noqa: E501
            # Mode
            self.mode = "cone"  # parallel, cone                ...
            self.filter = None
        if nVoxel is not None:
            # VARIABLE                                          DESCRIPTION                    UNITS
            # -------------------------------------------------------------------------------------
            self.DSD = 1536.0  # Distance Source Detector      (mm)
            self.DSO = 1000.0  # Distance Source Origin        (mm)
            # Detector parameters
            self.nDetector = np.array((nVoxel[1], nVoxel[2]))  # (V,U) number of pixels        (px)
            self.dDetector = np.array([0.8, 0.8])  # size of each pixel            (mm)
            self.sDetector = self.dDetector * self.nDetector  # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((nVoxel))  # number of voxels              (vx)
            self.sVoxel = np.array((256, 256, 256))  # total size of the image       (mm)
            self.dVoxel = self.sVoxel / self.nVoxel  # size of each voxel            (mm)
            # Offsets
            self.offOrigin = np.array((0, 0, 0))  # Offset of image from origin   (mm)
            self.offDetector = np.array((0, 0))  # Offset of Detector            (mm)
            self.rotDetector = np.array((0, 0, 0))
            # Auxiliary
            self.accuracy = 0.5  # Accuracy of FWD proj          (vx/sample)  # noqa: E501
            # Mode
            self.mode = "cone"  # parallel, cone
