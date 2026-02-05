from __future__ import division, absolute_import, print_function

# fix for DLL import issue on python 3.8 or higher on windows
# related to:
# https://github.com/CERN/TIGRE/issues/247
# https://github.com/CERN/TIGRE/issues/349
import os

_CUDA_PATH_WIN_DIR_HANDLE = None
def _add_cuda_path_on_windows():
    global _CUDA_PATH_WIN_DIR_HANDLE
    import sys
    if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
        _cuda_path = os.environ.get("CUDA_PATH")
        if _cuda_path:
            _bin = os.path.join(_cuda_path, "bin")
            if os.path.isdir(_bin):
                # a reference to the handle returned by os.add_dll_directory is necessary, 
                # otherwise it will be garbage collected and the directory will be removed 
                # from the search path again.
                _CUDA_PATH_WIN_DIR_HANDLE = os.add_dll_directory(_bin)

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

# always import Ax and _Ax_ext before all other ctypes extensions (checks for DLL import errors).
from .utilities.Ax import Ax, _Ax_ext, _ensure_Ax_ext_import
# _Ax_ext will be None if the import failed, 
# in that case we try to add the CUDA path on windows and import again, 
# this will raise the original import error with the full message if it still fails.
if _Ax_ext is None:
    from .utilities.Ax import _try_import_Ax_ext
    _add_cuda_path_on_windows()
    _try_import_Ax_ext()
_ensure_Ax_ext_import() # check if the import was successful, if not this will raise the original import error.

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
