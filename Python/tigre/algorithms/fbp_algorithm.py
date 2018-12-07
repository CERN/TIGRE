from __future__ import division
from __future__ import print_function
import os
import sys
import copy
from tigre.Atb import Atb
from tigre.utilities.filtering import filtering

currDir = os.path.dirname(os.path.realpath(__file__))
rootDir = os.path.abspath(os.path.join(currDir, '..'))
if rootDir not in sys.path:  # add parent dir to paths
    sys.path.append(rootDir)

def fbp(proj,geo,angles,filter=None,verbose=False):
    if geo.mode != 'parallel':
        raise ValueError("Only use FBP for parallel beam. Check geo.mode.")
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    proj_filt = filtering(proj,geox,angles,parker=False,verbose=verbose)

    res = Atb(proj_filt,geo,angles)*geo.DSO/geo.DSD

    return res






