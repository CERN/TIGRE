"""Object poses preserve hardware, transmission checks and reconstruction frames."""
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pytest
from conftest import poly_vertices
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from tigre.utilities.detector_geometry_builder import build_modular_surface
from tigre.utilities.geometry_io import save_geometry_yaml, load_geometry_file, surface_archive_arrays, load_geometry_archive
from tigre.utilities.object_motion import fixed_object_geometry, rotation_zyz, object_angles
from tigre.utilities.ct_physics import validate_transmission_geometry


def settings():
    return dict(layout="Piecewise linear", n_views=3, source_radius_mm=1000,
        source_z_mm=0, module_shape=(3, 4), pixel_spacing_mm=(3, 3),
        detector_radius_mm=536, column_gap=0, n_modules=1, angular_span_deg=0,
        transverse_centers_mm=[0], normal_offsets_mm=[0], panel_angles_deg=[0],
        z_offsets_mm=[0], fixed_hardware=True)


def volume():
    return SimpleNamespace(sVoxel=np.array([20., 40., 60.]), nVoxel=np.array([2, 4, 6]),
        offOrigin=np.array([1., 2., 3.]), object_angles_deg=np.array([[0, 0, 0], [30, 90, 60], [-20, 40, 90]]))


def test_fixed_hardware_and_inverse_rigid_transform():
    surface, geo = build_modular_surface(**settings()), volume()
    np.testing.assert_allclose(surface.source_positions_mm, np.tile(surface.source_positions_mm[0], (3, 1)))
    np.testing.assert_allclose(surface.pixel_centers_mm, np.repeat(surface.pixel_centers_mm[:1], 3, axis=0))
    transformed, fixed = fixed_object_geometry(surface, geo)
    for i, angles in enumerate(geo.object_angles_deg):
        rotation = rotation_zyz(angles)
        np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=1e-14)
        np.testing.assert_allclose(transformed.source_positions_mm[i] @ rotation.T + geo.offOrigin[::-1], surface.source_positions_mm[i], atol=1e-4)
        np.testing.assert_allclose(transformed.pixel_centers_mm[i] @ rotation.T + geo.offOrigin[::-1], surface.pixel_centers_mm[i], atol=1e-4)
    np.testing.assert_allclose(fixed.offOrigin, 0)
    assert validate_transmission_geometry(surface, geo)["crossing_samples"] > 0


@pytest.mark.parametrize("angles", [[], [[0, 1]], [[0, float('nan'), 0]], [[0, 0, 0]] * 2])
def test_invalid_angle_lists(angles):
    with pytest.raises(ValueError):
        object_angles(SimpleNamespace(object_angles_deg=angles), 3)




@pytest.mark.parametrize("angles", [[[30, 40, 50]], [[0, 0, 0], [30, 90, 60], [-20, 40, 90]], [0, 30, 90]])
def test_rotation_matlab_bridge_round_trip(tmp_path, angles):
    from tigre.utilities.geometry_matlab_bridge import to_mat, from_mat
    geo = volume()
    geo.object_angles_deg = np.asarray(angles)
    path, mat, output = (tmp_path / name for name in ("object.yaml", "object.mat", "output.yaml"))
    save_geometry_yaml(path, settings=settings(), geometry=geo)
    to_mat(path, mat)
    from_mat(mat, output)
    restored = load_geometry_file(output)
    np.testing.assert_allclose(object_angles(restored.geometry, 3), object_angles(geo, 3))
    assert restored.settings["fixed_hardware"]


def test_renderer_rotates_box_while_hardware_stays_fixed():
    from tigre.utilities.visualization.surface_renderer import SurfaceRenderer
    fig = plt.figure()
    try:
        renderer = SurfaceRenderer(fig.add_subplot(111, projection="3d"), build_modular_surface(**settings()), geo=volume())
        vertices = poly_vertices(renderer.volume).copy()
        source = np.asarray(renderer.source._offsets3d).copy()
        renderer.update(1)
        assert not np.allclose(vertices, poly_vertices(renderer.volume))
        np.testing.assert_allclose(renderer.source._offsets3d, source)
    finally:
        plt.close(fig)
