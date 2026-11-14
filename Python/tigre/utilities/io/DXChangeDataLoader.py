from __future__ import print_function
from __future__ import with_statement

import os
import math
import numpy as np
from tqdm import tqdm
import h5py

from configparser import ConfigParser
from tigre.utilities.geometry import geometry

def DXChangeDataLoader(filename):
    """Load a DXchange (HDF5) dataset, the format used at many synchrotrons.

    Flat- and dark-field corrects the projections, converts them to line
    integrals and builds a parallel-beam geometry with 1 mm voxels and pixels
    (rescale ``geo`` to your pixel size if you need physical units).

    Parameters
    ----------
    filename : str
        Path of the HDF5 file, with ``/exchange/data``, ``data_white``,
        ``data_dark`` and ``theta`` (in degrees).

    Returns
    -------
    proj : numpy.ndarray
        Projections (line integrals), float32, shape ``(N, v, u)``.
    geo : tigre.utilities.geometry.Geometry
        A parallel-beam geometry matching the data.
    angles : numpy.ndarray
        Projection angles in radians.
    """

    with h5py.File(filename, "r") as f:
        group = f['/exchange/']
        proj = group['data'][()]
        white = group['data_white'][()]
        dark = group['data_dark'][()]
        angles = group['theta'][()]*np.pi/180
    f.close()

    white = np.float32(white)
    dark = np.float32(dark)
    proj = np.float32(proj)

    white[white<0] = 0
    dark[dark<0] = 0
    dark = np.mean(dark,axis=0)
    white = np.mean(white,axis=0)
    proj[proj<=0]=0.001
    proj = (proj-dark)/(white-dark)
    proj = -np.log(proj/proj.max())
    geo = geometry(mode='parallel',nVoxel=np.array([proj.shape[1], proj.shape[2], proj.shape[2]]))

    return np.float32(proj), geo, angles