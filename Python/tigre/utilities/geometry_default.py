from __future__ import division

import numpy as np
from tigre.utilities.geometry import Geometry


class ConeGeometryDefault(Geometry):
    """A ready-made cone-beam geometry, for tests and examples.

    Available as ``tigre.geometry_default``. Source to detector 1536 mm,
    source to origin 1000 mm, a 256 mm cube image centred on the origin.

    Parameters
    ----------
    high_resolution : bool, optional
        True (default): 512³ voxels and a 512 x 512 detector of 0.8 mm pixels.
        False: 64³ voxels and a 128 x 128 detector of 3.2 mm pixels.
    nVoxel : numpy.ndarray, optional
        Image size ``(z, y, x)``; the detector then has ``(nVoxel[1],
        nVoxel[2])`` pixels of 0.8 mm. Overrides ``high_resolution``.

    Examples
    --------
    >>> import tigre
    >>> geo = tigre.geometry_default(high_resolution=False)
    >>> geo.nVoxel
    array([64, 64, 64])
    """

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

class FanGeometryDefault(Geometry):
    """A ready-made 2D fan-beam geometry: a one-slice cone-beam geometry.

    Available as ``tigre.fan_geometry_default``.

    Parameters
    ----------
    high_resolution : bool, optional
        True (default): a 512 x 512 image and 512 detector pixels. False: a
        64 x 64 image and 128 detector pixels.
    nVoxel : numpy.ndarray, optional
        Image size ``(1, y, x)``; must have one slice.

    Raises
    ------
    ValueError
        If ``nVoxel`` has more than one slice.
    """

    def __init__(self, high_resolution=True, nVoxel=None):
        if nVoxel is not None and nVoxel[0] !=1:
            raise ValueError("Fan beam geometry only supports 2D images")
    
        Geometry.__init__(self)
        if high_resolution:
            # VARIABLE                                          DESCRIPTION                    UNITS
            # -------------------------------------------------------------------------------------
            self.DSD = 1536.0  # Distance Source Detector      (mm)
            self.DSO = 1000.0  # Distance Source Origin        (mm)
            # Detector parameters
            self.nDetector = np.array((1, 512))  # number of pixels              (px)
            self.dDetector = np.array((0.8, 0.8))  # size of each pixel            (mm)
            self.sDetector = self.nDetector * self.dDetector  # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((1, 512, 512))  # number of voxels              (vx)
            self.sVoxel = np.array((0.5, 256, 256))  # total size of the image       (mm)
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
            self.nDetector = np.array((1, 128))  # (V,U) number of pixels        (px)
            self.dDetector = np.array((0.8, 0.8)) * 4  # size of each pixel            (mm)
            self.sDetector = self.nDetector * self.dDetector  # total size of the detector    (mm)
            # Image parameters
            self.nVoxel = np.array((1, 64, 64))  # number of voxels              (vx)
            self.sVoxel = np.array((0.5, 256, 256))  # total size of the image       (mm)
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
            
            # Image parameters
            self.nVoxel = np.array((nVoxel))  # number of voxels              (vx)
            self.sVoxel = np.array((256, 256, 256))  # total size of the image       (mm)
            self.dVoxel = self.sVoxel / self.nVoxel  # size of each voxel            (mm)
            # Detector parameters
            self.nDetector = np.array((nVoxel[0], nVoxel[1]))  # (V,U) number of pixels        (px)
            self.sDetector = self.sVoxel[1:]*self.DSD/self.DSO
            self.dDetector = self.sDetector / self.nDetector  # size of each pixel            (mm)
            # Offsets
            self.offOrigin = np.array((0, 0, 0))  # Offset of image from origin   (mm)
            self.offDetector = np.array((0, 0))  # Offset of Detector            (mm)
            self.rotDetector = np.array((0, 0, 0))
            # Auxiliary
            self.accuracy = 0.5  # Accuracy of FWD proj          (vx/sample)  # noqa: E501
            # Mode
            self.mode = "cone"  # parallel, cone
 