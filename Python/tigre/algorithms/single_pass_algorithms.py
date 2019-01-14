from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import sys
import numpy as np
import copy
from tigre.utilities.Atb import Atb
from tigre.utilities.filtering import filtering

# TODO: this is quite nasty; it would be nice to reorganise file structure later so top level folder is always in path
currDir = os.path.dirname(os.path.realpath(__file__))
rootDir = os.path.abspath(os.path.join(currDir, '..'))
if rootDir not in sys.path:  # add parent dir to paths
    sys.path.append(rootDir)
def FDK(proj, geo, angles, filter=None,verbose=False):
    ('\n'
     'FDK solves Cone Beam CT image reconstruction'
     '\n'
     'Parameters \n'
     '-------------------------------------------------------------------\n'
     '\n'
     'proj:         Data input in the form of 3d, np.array(dtype=float32)\n'
     '\n'
     'geo:          Geometry of detector and image (see examples/Demo code)\n'
     '\n'
     'alpha:        1d array of angles corresponding to image data, np.array(dtype=float32)\n'
     '\n'
     'filter:       default="ram_lak" \n'
     '              opts: \n'
     '                  "shep_logan"'
     '                  "cosine"   '
     '                  "hamming" '
     '                  "hann" '
     'Examples \n'
     '---------------------------------------------------------------------------\n'
     'See Demos/Example code for further instructions.\n'
     '---------------------------------------------------------------------------'
     '\n'
     """This file is part of the TIGRE Toolbox

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
                           PYTHON : Reuben Lindroos,Sam Loescher, """)
    geo = copy.deepcopy(geo)
    geo.check_geo(angles)
    if filter is not None:
        geo.filter = filter
    # Weight
    proj_filt = np.zeros(proj.shape, dtype=np.float32)
    for ii in range(angles.shape[0]):
        xv = np.arange((-geo.nDetector[1] / 2) + 0.5, 1 + (geo.nDetector[1] / 2) - 0.5) * geo.dDetector[1]
        yv = np.arange((-geo.nDetector[0] / 2) + 0.5, 1 + (geo.nDetector[0] / 2) - 0.5) * geo.dDetector[0]
        (xx, yy) = np.meshgrid(xv, yv)

        w = geo.DSD[0] / np.sqrt((geo.DSD[0] ** 2 + xx ** 2 + yy ** 2))
        proj_filt[ii] = proj[ii] * w

    proj_filt = filtering(proj_filt, geo, angles, parker=False,verbose=verbose)
    # m = {
    #     'py_projfilt': proj_filt,
    #
    # }
    # scipy.io.savemat('Tests/Filter_data', m)
    res = Atb(proj_filt, geo, geo.angles, 'FDK')
    # res = 0
    # res = Atb(proj,geo,angles,'FDK')
    return res


fdk = FDK

def fbp(proj,geo,angles,filter=None,verbose=False):
    if geo.mode != 'parallel':
        raise ValueError("Only use FBP for parallel beam. Check geo.mode.")
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    proj_filt = filtering(proj,geox,angles,parker=False,verbose=verbose)

    res = Atb(proj_filt,geo,angles)*geo.DSO/geo.DSD

    return res