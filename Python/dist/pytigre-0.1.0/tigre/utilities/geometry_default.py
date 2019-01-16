from __future__ import division
import numpy as np
from tigre.utilities.geometry import Geometry

class TIGREParameters(Geometry):

    def __init__(self, high_quality=True):
        if high_quality:
            # VARIABLE                                          DESCRIPTION                    UNITS
            # -------------------------------------------------------------------------------------
            self.DSD = 1536                                     # Distance Source Detector      (mm)
            self.DSO = 1000                                     # Distance Source Origin        (mm)
            # Detector parameters
            self.nDetector = np.array((512, 512))               # number of pixels              (px)
            self.dDetector = np.array((0.8, 0.8))               # size of each pixel            (mm)
            self.sDetector = self.nDetector * self.dDetector    # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((256, 256, 256))             # number of voxels              (vx)
            self.sVoxel = np.array((256, 256, 256))             # total size of the image       (mm)
            self.dVoxel = self.sVoxel/self.nVoxel               # size of each voxel            (mm)
            # Offsets
            self.offOrigin = np.array((0, 0, 0))                # Offset of image from origin   (mm)
            self.offDetector = np.array((0, 0))                 # Offset of Detector            (mm)

            # Auxiliary
            self.accuracy = 0.5                                 # Accuracy of FWD proj          (vx/sample)
            # Mode
            self.mode = 'cone'                                  # parallel, cone                ...
            self.filter = None
        else:
            # VARIABLE                                          DESCRIPTION                    UNITS
            # -------------------------------------------------------------------------------------
            self.DSD = 1536                                     # Distance Source Detector      (mm)
            self.DSO = 1000                                     # Distance Source Origin        (mm)
            # Detector parameters
            self.nDetector = np.array((128, 128))                # (V,U) number of pixels        (px)
            self.dDetector = np.array((0.8, 0.8))*4             # size of each pixel            (mm)
            self.sDetector = self.nDetector * self.dDetector    # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((64, 63 , 62))               # number of voxels              (vx)
            self.sVoxel = np.array((256, 256, 256))             # total size of the image       (mm)
            self.dVoxel = self.sVoxel / self.nVoxel             # size of each voxel            (mm)
            # Offsets
            self.offOrigin = np.array((0, 0, 0))                # Offset of image from origin   (mm)
            self.offDetector = np.array((0, 0))                 # Offset of Detector            (mm)
            self.rotDetector = np.array((0,0,0))
            # Auxiliary
            self.accuracy = 0.5                                 # Accuracy of FWD proj          (vx/sample)
            # Mode
            self.mode='cone'                                    # parallel, cone                ...
            self.filter=None
        Geometry.__init__(self)
    def __str__(self):
        parameters = []
        parameters.append("TIGRE parameters")
        parameters.append("-----")
        parameters.append("Distance from source to detector = " + str(self.DSD) + "mm")
        parameters.append("Distance from source to origin = " + str(self.DSO) + "mm")

        parameters.append("-----")
        parameters.append("Detector parameters")
        parameters.append("Number of pixels = " + str(self.nDetector))
        parameters.append("Size of each pixel = " + str(self.dDetector) + "mm")
        parameters.append("Total size of the detector = " + str(self.sDetector) + "mm")

        parameters.append("-----")
        parameters.append("Image parameters")
        parameters.append("Number of voxels = " + str(self.nVoxel))
        parameters.append("Total size of the image = " + str(self.sVoxel) + "mm")
        parameters.append("Size of each voxel = " + str(self.dVoxel) + "mm")

        parameters.append("-----")
        parameters.append("Offset correction parameters")
        parameters.append("Offset of image from origin = " + str(self.offOrigin) + "mm")
        parameters.append("Offset of detector = " + str(self.offDetector) + "mm")

        parameters.append("-----")
        parameters.append("Auxillary parameters")
        parameters.append("Accuracy of forward projection = " + str(self.accuracy))

        return '\n'.join(parameters)
