#%% Demo 23: Curved detector
#
# This demo higlights the curved detector processing.
#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# This file is part of the TIGRE Toolbox
#
# Copyright (c) 2015, University of Bath and
#                     CERN-European Organization for Nuclear Research
#                     All rights reserved.
#
# License:            Open Source under BSD.
#                     See the full license at
#                     https://github.com/CERN/TIGRE/blob/master/LICENSE
#
# Contact:            tigre.toolbox@gmail.com
# Codes:              https://github.com/CERN/TIGRE/
# Coded by:           Ander Biguri
# --------------------------------------------------------------------------
#%%Initialize
import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs

#%% Geometry
geo = tigre.geometry_default(high_resolution=False)

# If you have a curved detector, currently the only way to process it is to "flatten" the detector. 
from tigre.utilities.curved_detector import flatten_detector

proj= ...
# You can change the oversample value to create denser (and thus more accurate) flat projections, in exchange for computational time at recon
proj = flatten_detector(proj,geo,oversample=1)


# Now just use as if proj came from flat panel detectors. 