from __future__ import division
import tigre
import numpy as np

def geo(tomoshape, **kwargs):
    """
    Define TIGRE geo according to shape of object being reconstructed.
    Parameters:
    ----------
    tomoshape: np.array()
    shape of obj being reconstructed as array

    Extra options
    ----------
    DSD : int                                   # Distance Source Detector      (mm)
    DSO : int                                   # Distance Source Origin        (mm)

    # Detector parameters

    nDetector : np.array()                      # number of pixels              (px)
    dDetector : np.array((0.8, 0.8))*4          # size of each pixel            (mm)
    sDetector : nDetector * dDetector           # total size of the detector    (mm)

    # Image parameters

    nVoxel : np.array()                         # number of voxels              (vx)
    sVoxel : np.array()                         # total size of the image       (mm)
    dVoxel : sVoxel / nVoxel                    # size of each voxel            (mm)

    # Offsets

    offOrigin : np.array()                      # Offset of image from origin   (mm)
    offDetector : np.array()                    # Offset of Detector            (mm)

    # Auxiliary

    accuracy : float                            # Accuracy of FWD proj          (vx/sample)

    # Mode
    mode : str                                   # parallel, cone                ...
    filter : str                                 # FDK filter to user (see help(FDK))

    Example
    ---------
    >>> import numpy as np
    >>> import tigre
    >>> obj = np.zeros((64,128,256))
    >>> geo = tigre.geo(obj,geo={'filter':'ram_lak','nVoxel':np.array((64,64,64))})
    >>> print(geo.nVoxel, geo.filter)
    # returns: (array([64, 64, 64]), 'ram_lak')

    See example code for futher help.

     ---------------------------------------------------------------------------
        This file is part of the TIGRE Toolbox

        Copyright (c) 2015, University of Bath and
                            CERN-European Organization for Nuclear Research
                            All rights reserved.

        License:            Open Source under BSD.
                            See the full license at
                            https://github.com/CERN/TIGRE/license.txt

        Contact:            tigre.toolbox@gmail.com
        Codes:              https://github.com/CERN/TIGRE/
      --------------------------------------------------------------------------
        Coded by:          MATLAB (original code): Ander Biguri
                           PYTHON : Sam Loescher, Reuben Lindroos"""

    from tigre.geometry import TIGREParameters
    geox = TIGREParameters(high_quality=False)
    default_opts = {'DSD': 1536,
                    'DSO': 1000,
                    'nDetector': np.array((128, 128)),
                    'dDetector': np.array((0.8, 0.8)) * 4,
                    'nVoxel': tomoshape,
                    'sVoxel': np.array((256, 256, 256)),
                    'offOrigin': np.array((0, 0, 0)),
                    'offDetector': np.array((0, 0)),
                    'accuracy': 0.5,
                    'mode': None,
                    'filter': None
                    }
    geo = default_opts

    if 'geo' in kwargs:
        geo.update(kwargs['geo'])

    geox.DSD = geo['DSD']
    geox.DSO = geo['DSO']
    # Detector parameters
    geox.nDetector = geo['nDetector']
    geox.dDetector = geo['dDetector']
    geox.sDetector = geox.nDetector * geox.dDetector
    # Image parameters
    geox.nVoxel = geo['nVoxel']
    geox.sVoxel = geo['sVoxel']
    geox.dVoxel = geox.sVoxel / geox.nVoxel
    # Offsets
    geox.offOrigin = geo['offOrigin']
    geox.offDetector = geo['offDetector']

    # Auxiliary
    geox.accuracy = geo['accuracy']
    # Mode
    geox.mode = geo['mode']
    geox.filter = geo['filter']
    return geox