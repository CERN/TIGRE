"""
Curved detectors
================

Resample curved-detector data onto a flat detector, then reconstruct.

Source: `Python/demos/d24_CurvedDetector.py <https://github.com/CERN/TIGRE/blob/master/Python/demos/d24_CurvedDetector.py>`__ · `MATLAB/Demos/d024_CurvedDetector.m <https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d024_CurvedDetector.m>`__

.. dropdown:: MATLAB version of this demo
   :icon: code

   .. literalinclude:: /../MATLAB/Demos/d024_CurvedDetector.m
      :language: matlab

This demo is not run when the website is built.
"""

# %%
# Initialize
# ----------

import tigre
import numpy as np
from tigre.utilities import sample_loader
from tigre.utilities import CTnoise
import tigre.algorithms as algs

# %%
# Geometry
# --------

geo = tigre.geometry_default(high_resolution=False)

# If you have a curved detector, currently the only way to process it is to "flatten" the detector.
from tigre.utilities.curved_detector import flatten_detector

proj= ...
# You can change the oversample value to create denser (and thus more accurate) flat projections, in exchange for computational time at recon
proj = flatten_detector(proj,geo,oversample=1)

# Now just use as if proj came from flat panel detectors.
