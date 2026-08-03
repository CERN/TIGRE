"""Regression coverage for Varian geometry conversion."""

from contextlib import contextmanager
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import numpy as np


@contextmanager
def _varian_loader_module():
    """Load the Varian I/O module without requiring TIGRE's CUDA extensions."""
    package_root = Path(__file__).resolve().parents[1] / "tigre"
    package_paths = {
        "tigre": package_root,
        "tigre.utilities": package_root / "utilities",
        "tigre.utilities.io": package_root / "utilities" / "io",
        "tigre.utilities.io.varian": package_root / "utilities" / "io" / "varian",
    }
    saved_modules = {
        name: module
        for name, module in sys.modules.items()
        if name == "tigre" or name.startswith("tigre.")
    }

    for name in saved_modules:
        del sys.modules[name]
    for name, path in package_paths.items():
        package = ModuleType(name)
        package.__path__ = [str(path)]
        sys.modules[name] = package

    try:
        yield importlib.import_module("tigre.utilities.io.varian.varian_io")
    finally:
        for name in list(sys.modules):
            if name == "tigre" or name.startswith("tigre."):
                del sys.modules[name]
        sys.modules.update(saved_modules)


def test_fallback_geometry_sizes_lateral_offset():
    scan_params = SimpleNamespace(
        SID=1500,
        SAD=1000,
        imager_size=[400, 300],
        imager_res=[0.5, 0.4],
        imager_lat=20,
    )
    zero_offset_scan_params = SimpleNamespace(**{**scan_params.__dict__, "imager_lat": 0})

    with _varian_loader_module() as varian_io:
        offset_geometry = varian_io.read_varian_geometry(scan_params, recon_params=None)
        zero_offset_geometry = varian_io.read_varian_geometry(
            zero_offset_scan_params, recon_params=None
        )

    np.testing.assert_array_equal(offset_geometry.offDetector, [0, -20])
    np.testing.assert_array_equal(offset_geometry.nVoxel, [400, 340, 340])
    np.testing.assert_array_equal(zero_offset_geometry.nVoxel, [400, 300, 300])
