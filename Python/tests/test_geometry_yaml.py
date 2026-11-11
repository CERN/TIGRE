"""YAML round trips, exact companion surfaces, animation and GUI loading."""
import numpy as np
import pytest
from conftest import poly_vertices
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import yaml
from pathlib import Path
from unittest.mock import patch
from tigre.utilities.geometry_io import save_geometry_yaml, load_geometry_file
from tigre.utilities.detector_geometry_builder import build_modular_surface
from tigre.utilities.detector_geometry import DetectorSurfaceGeometry



def settings(layout="Piecewise linear", motion="Rotating source + detector"):
    return dict(layout=layout, n_views=3, source_radius_mm=1000, source_z_mm=5,
        module_shape=(3, 4), pixel_spacing_mm=(2, 3), detector_radius_mm=536,
        column_gap=1, n_modules=3, angular_span_deg=40,
        transverse_centers_mm=[-20, 0, 20], normal_offsets_mm=[0, 5, 0],
        panel_angles_deg=[10, 0, -10], z_offsets_mm=[0, 3, 0],
        panel_x_angles_deg=[5, 0, -5], panel_y_angles_deg=[0, 10, -10],
        source_motion=motion, source_start_mm=[1000, -100, 0],
        source_end_mm=[1000, 100, 20], orbit_start_deg=-30, scan_span_deg=60,
        cor_mm=7 if motion == "Rotating source + detector" else 0)


def test_moving_volume_yaml_and_renderer(tmp_path):
    from types import SimpleNamespace
    from tigre.utilities.visualization.surface_renderer import SurfaceRenderer
    offsets = np.array([[-100, 2, 3], [0, 2, 3], [100, 2, 3]])
    geo = SimpleNamespace(sVoxel=np.array([20, 30, 40]),
                          nVoxel=np.array([10, 15, 20]), offOrigin=offsets)
    path = tmp_path / "moving.yaml"
    save_geometry_yaml(path, settings=settings(), geometry=geo)
    loaded = load_geometry_file(path)
    np.testing.assert_array_equal(loaded.geometry.offOrigin, offsets)
    np.testing.assert_array_equal(loaded.geometry.nVoxel, geo.nVoxel)
    fig = plt.figure()
    try:
        ax = fig.add_subplot(111, projection="3d")
        renderer = SurfaceRenderer(ax, loaded.surface, geo=loaded.geometry)
        initial = poly_vertices(renderer.volume).copy()
        renderer.update(2)
        np.testing.assert_allclose(poly_vertices(renderer.volume)[2] - initial[2], 200)
        np.testing.assert_allclose(poly_vertices(renderer.volume)[:2], initial[:2])
        assert ax.get_zlim3d()[0] < -110
        assert ax.get_zlim3d()[1] > 110
    finally:
        plt.close(fig)


@pytest.mark.parametrize("layout", ["Piecewise linear", "Cylindrical"])
@pytest.mark.parametrize("motion", ["Rotating source + detector",
    "Fixed detector, moving source", "Fixed detector, arc source"])
def test_editable_yaml_round_trip(layout, motion, tmp_path):
    parameters = settings(layout, motion)
    expected = build_modular_surface(**parameters)
    path = tmp_path / "geometry.yaml"
    save_geometry_yaml(path, settings=parameters)
    assert "!!python" not in path.read_text()
    assert not list(tmp_path.glob("*.npz"))
    loaded = load_geometry_file(path)
    assert loaded.settings["layout"] == layout
    assert loaded.settings["source_motion"] == motion
    np.testing.assert_array_equal(loaded.surface.pixel_centers_mm, expected.pixel_centers_mm)
    np.testing.assert_array_equal(loaded.surface.source_positions_mm, expected.source_positions_mm)
    np.testing.assert_array_equal(loaded.surface.valid_mask, expected.valid_mask)
    assert loaded.surface._module_slices == expected._module_slices


def test_exact_surface_yaml_companion_and_animation(tmp_path):
    from types import SimpleNamespace
    original = build_modular_surface(**settings())
    mask = original.valid_mask.copy()
    mask[1, 0, 0] = False
    surface = DetectorSurfaceGeometry(original.source_positions_mm,
        original.pixel_centers_mm, mask, module_slices=original._module_slices)
    volume = SimpleNamespace(sVoxel=np.array([20, 30, 40]), offOrigin=np.array([1, 2, 3]))
    path = tmp_path / "geometry.yml"
    save_geometry_yaml(path, surface=surface, geometry=volume)
    data = yaml.safe_load(path.read_text())
    assert data["kind"] == "surface"
    assert not Path(data["surface_file"]).is_absolute()
    loaded = load_geometry_file(path)
    assert loaded.settings is None
    np.testing.assert_array_equal(loaded.surface.pixel_centers_mm, surface.pixel_centers_mm)
    np.testing.assert_array_equal(loaded.surface.valid_mask, surface.valid_mask)
    assert loaded.surface._module_slices == surface._module_slices
    np.testing.assert_array_equal(loaded.geometry.sVoxel, volume.sVoxel)





@pytest.mark.parametrize("change", [
    {"version": 99}, {"units": {"length": "cm", "angle": "degree"}},
    {"rotation_convention": "unknown"}, {"kind": "unknown"}, {"typo": 1},
    {"builder": {"unknown": 1}},
])
def test_yaml_rejects_incompatible_or_malformed_settings(change, tmp_path):
    path = tmp_path / "geometry.yaml"
    save_geometry_yaml(path, settings=settings())
    data = yaml.safe_load(path.read_text())
    data.update(change)
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError):
        load_geometry_file(path)


def test_yaml_does_not_construct_python_objects(tmp_path):
    path = tmp_path / "unsafe.yaml"
    path.write_text("!!python/object:builtins.object {}")
    with pytest.raises(ValueError, match="Invalid geometry YAML"):
        load_geometry_file(path)







def test_save_rejects_volume_offsets_with_wrong_view_count(tmp_path):
    from types import SimpleNamespace
    volume = SimpleNamespace(sVoxel=np.array([20, 30, 40]), offOrigin=np.zeros((2, 3)))
    with pytest.raises(ValueError, match="per view"):
        save_geometry_yaml(tmp_path / "wrong.yaml", settings=settings(), geometry=volume)


def test_failed_atomic_save_keeps_existing_yaml(tmp_path):
    path = tmp_path / "geometry.yaml"
    save_geometry_yaml(path, settings=settings())
    original = path.read_bytes()
    with patch("pathlib.Path.replace", side_effect=OSError("write blocked")):
        with pytest.raises(OSError):
            save_geometry_yaml(path, settings=settings(layout="Cylindrical"))
    assert path.read_bytes() == original
    assert not list(tmp_path.glob("*.tmp"))
