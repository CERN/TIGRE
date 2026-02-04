from __future__ import division, absolute_import, print_function

# fix for DLL import issue on python 3.8 or higher on windows
# related to:
# https://github.com/CERN/TIGRE/issues/247
# https://github.com/CERN/TIGRE/issues/349
# https://github.com/CERN/TIGRE/issues/716
import os
import sys

_DLL_DIR_HANDLES = []
if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
    _cuda_path = os.environ.get("CUDA_PATH")
    if _cuda_path:
        _bin = os.path.join(_cuda_path, "bin")
        if os.path.isdir(_bin):
            _DLL_DIR_HANDLES.append(os.add_dll_directory(_bin))

# if hasattr(os, "add_dll_directory"):
#     # Add all the DLL directories manually
#     # see:
#     # https://docs.python.org/3.8/whatsnew/3.8.html#bpo-36085-whatsnew
#     # https://stackoverflow.com/a/60803169/19344391
# dll_directory = os.path.dirname(__file__)
# os.add_dll_directory(dll_directory)

#     # The user must install the CUDA Toolkit
# cuda_bin = os.path.join(os.environ["CUDA_PATH"], "bin")
# os.add_dll_directory(cuda_bin)

from .utilities.geometry import geometry
from .utilities.geometry_default import ConeGeometryDefault as geometry_default
from .utilities.geometry_default import FanGeometryDefault as fan_geometry_default
from .utilities.Ax import Ax
from .utilities.Atb import Atb
from .utilities.visualization.plotproj import plotproj, plotProj, plotSinogram
from .utilities.visualization.plotimg import plotimg, plotImg
from .utilities.visualization.plot_geometry import plot_geometry
from .utilities.visualization.plot_angles import plot_angles
from .utilities.CTnoise import add
from .utilities.common_geometry import (
    staticDetectorGeo,
    staticDetLinearSourceGeo,
    ArbitrarySourceDetMoveGeo,
)

from . import algorithms
