
Demo 01: Describing your geometry
=================================

To see a demo of what the geometry paramterers should look like, do as
follows:

.. code:: ipython2

    import tigre
    geo = tigre.geometry_default(high_quality = False)
    print(geo)



.. parsed-literal::

    TIGRE parameters
    -----
    Geometry parameters
    Distance from source to detector (DSD) = 1536 mm
    Distance from source to origin (DSO)= 1000 mm
    -----
    Detector parameters
    Number of pixels (nDetector) = [128 128]
    Size of each pixel (dDetector) = [3.2 3.2] mm
    Total size of the detector (sDetector) = [409.6 409.6] mm
    -----
    Image parameters
    Number of voxels (nVoxel) = [64 64 64]
    Total size of the image (sVoxel) = [256 256 256] mm
    Size of each voxel (dVoxel) = [4. 4. 4.] mm
    -----
    Offset correction parameters
    Offset of image from origin (offOrigin) = [0 0 0] mm
    Offset of detector (offDetector) = [0 0] mm
    -----
    Auxillary parameters
    Samples per pixel of forward projection (accuracy) = 0.5
    -----
    Rotation of the Detector (rotDetector) = [0 0 0] rad


We recommend using the template below and defining you’re class as such:

.. code:: ipython2

    from __future__ import division
    import numpy as np
    
    
    class Geometry:
    
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
            else:
                # VARIABLE                                          DESCRIPTION                    UNITS
                # -------------------------------------------------------------------------------------
                self.DSD = 1536                                     # Distance Source Detector      (mm)
                self.DSO = 1000                                     # Distance Source Origin        (mm)
                # Detector parameters
                self.nDetector = np.array((128, 128))             # number of pixels              (px)
                self.dDetector = np.array((0.8, 0.8))*4             # size of each pixel            (mm)
                self.sDetector = self.nDetector * self.dDetector    # total size of the detector    (mm)
                # Image parameters
                self.nVoxel = np.array((64, 64 , 64))           # number of voxels              (vx)
                self.sVoxel = np.array((256, 256, 256))             # total size of the image       (mm)
                self.dVoxel = self.sVoxel / self.nVoxel             # size of each voxel            (mm)
                # Offsets
                self.offOrigin = np.array((0, 0, 0))                # Offset of image from origin   (mm)
                self.offDetector = np.array((0, 0))                 # Offset of Detector            (mm)
    
                # Auxiliary
                self.accuracy = 0.5                                 # Accuracy of FWD proj          (vx/sample)
                # Mode
                self.mode=None                                      # parallel, cone                ...
                self.filter=None
