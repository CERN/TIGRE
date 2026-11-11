"""Core geometry files must not depend on the optional desktop application."""

import subprocess
import sys


def test_geometry_round_trips_without_gui_imports(tmp_path):
    script = """
import importlib.abc
import sys
from pathlib import Path

class ForbidGui(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "tkinter" or fullname.startswith("tkinter.") or (
            fullname == "tigre_gui" or fullname.startswith("tigre_gui.")
        ) or (
            fullname == "tigre_geometry_design_tool"
            or fullname.startswith("tigre_geometry_design_tool.")
        ) or fullname == "tigre.utilities.detector_geometry_gui":
            raise AssertionError("Unexpected GUI import: " + fullname)
        return None

sys.meta_path.insert(0, ForbidGui())
import numpy as np
from tigre.utilities.detector_geometry_builder import build_modular_surface
from tigre.utilities.geometry_io import load_geometry_file, save_geometry_yaml

settings = dict(layout="Piecewise linear", n_views=3, source_radius_mm=1000,
    source_z_mm=0, module_shape=(2, 3), pixel_spacing_mm=(1, 1),
    detector_radius_mm=500, column_gap=0, n_modules=1, angular_span_deg=0,
    transverse_centers_mm=[0], normal_offsets_mm=[0], panel_angles_deg=[0],
    z_offsets_mm=[0])
surface = build_modular_surface(**settings)
folder = Path(sys.argv[1])
recipe = folder / "recipe.yaml"
exact = folder / "exact.yaml"
save_geometry_yaml(recipe, settings=settings)
save_geometry_yaml(exact, surface=surface)
for path in (recipe, exact):
    loaded = load_geometry_file(path).surface
    np.testing.assert_array_equal(loaded.pixel_centers_mm, surface.pixel_centers_mm)
    np.testing.assert_array_equal(loaded.source_positions_mm, surface.source_positions_mm)
    np.testing.assert_array_equal(loaded.valid_mask, surface.valid_mask)
"""
    subprocess.run([sys.executable, "-c", script, str(tmp_path)], check=True)
