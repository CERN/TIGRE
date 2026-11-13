"""Acquisition cadence, multi-turn orbits and synchronized volume motion."""
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from tigre.utilities.detector_geometry_builder import acquisition_timing, build_modular_surface
from tigre.utilities.geometry_io import save_geometry_yaml, load_geometry_file, surface_archive_arrays, load_geometry_archive


def timing(**changes):
    return dict(readout_period_ms=100, firing_rate_hz=4, rotation_period_s=1,
                table_speed_mm_s=-10, playback_rate=2, **changes)


def parameters():
    return dict(layout="Piecewise linear", n_views=9, source_radius_mm=1000,
                source_z_mm=0, module_shape=(2, 3), pixel_spacing_mm=(1, 1),
                detector_radius_mm=500, column_gap=0, n_modules=1,
                angular_span_deg=0, transverse_centers_mm=[0], normal_offsets_mm=[0],
                panel_angles_deg=[0], z_offsets_mm=[0], acquisition=timing())






@pytest.mark.parametrize("field,value", [
    ("readout_period_ms", 0), ("readout_period_ms", 300),
    ("firing_rate_hz", -1), ("rotation_period_s", 0),
    ("table_speed_mm_s", np.nan), ("playback_rate", 0),
])
def test_invalid_timing_is_rejected(field, value):
    values = timing()
    values[field] = value
    with pytest.raises(ValueError):
        acquisition_timing(9, values)
