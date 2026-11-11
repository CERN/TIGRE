"""Numeric loading and Tk inspection regression tests (no GPU required)."""
import numpy as np
import pytest
from unittest.mock import patch
from tigre.utilities.detector_geometry import DetectorSurfaceGeometry, surface_from_geometry
from tigre.utilities.geometry_io import load_geometry_archive, surface_archive_arrays


def sample_surface():
    p = np.zeros((3, 3, 6, 3))
    p[..., 0] = -100
    p[..., 1] = np.arange(6)
    p[..., 2] = np.arange(3)[None, :, None]
    p[1] += [0, 5, 3]
    mask = np.ones(p.shape[:3], dtype=bool)
    mask[1, 0, 0] = False
    return DetectorSurfaceGeometry([[100, 0, 0], [0, 100, 20], [-100, 0, 0]], p, mask,
        module_slices=[(slice(0, 3), slice(0, 3)), (slice(0, 3), slice(3, 6))])


def test_surface_archive_preserves_views_masks_modules_and_volume(tmp_path):
    from types import SimpleNamespace
    surface = sample_surface()
    volume = SimpleNamespace(sVoxel=np.array([20, 30, 40]), offOrigin=np.array([1, 2, 3]),
                             nVoxel=np.array([10, 15, 8]))
    path = tmp_path / "surface.npz"
    np.savez_compressed(path, **surface_archive_arrays(surface, volume))
    geo, loaded = load_geometry_archive(path)
    for name in ("source_positions_mm", "pixel_centers_mm", "valid_mask"):
        np.testing.assert_array_equal(getattr(loaded, name), getattr(surface, name))
    assert loaded._module_slices == surface._module_slices
    np.testing.assert_array_equal(geo.sVoxel, volume.sVoxel)
    np.testing.assert_array_equal(geo.offOrigin, volume.offOrigin)
    np.testing.assert_array_equal(geo.nVoxel, volume.nVoxel)
    np.testing.assert_allclose(geo.dVoxel, [2, 2, 5])


@pytest.mark.parametrize("counts", [[0, 2, 3], [1.5, 2, 3], [1, 2], [1, np.inf, 3]])
def test_volume_rejects_invalid_voxel_counts(counts):
    from tigre.utilities.geometry_io import _volume
    with pytest.raises(ValueError):
        _volume(dict(sVoxel=[20, 30, 40], nVoxel=counts))


def test_standard_archive_converts_per_view_tilt_and_cor(tmp_path):
    import tigre
    geo = tigre.geometry_default(high_resolution=False)
    geo.nDetector = np.array([3, 4])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.rotDetector = np.deg2rad([[0, 0, 0], [10, 15, 20], [-10, -15, -20]])
    geo.COR = np.array([0, 10, -10])
    angles = np.array([0, 0.5, 1.2])
    arrays = {name: getattr(geo, name) for name in (
        "DSO", "DSD", "nDetector", "dDetector", "rotDetector", "COR")}
    path = tmp_path / "standard.npz"
    np.savez(path, angles=angles, **arrays)
    _, actual = load_geometry_archive(path)
    expected = surface_from_geometry(geo, angles)
    np.testing.assert_allclose(actual.pixel_centers_mm, expected.pixel_centers_mm)
    np.testing.assert_allclose(actual.source_positions_mm, expected.source_positions_mm)


def test_old_surface_archive_and_bad_module_metadata(tmp_path):
    surface = sample_surface()
    arrays = surface_archive_arrays(surface)
    arrays.pop("module_slices")
    path = tmp_path / "old.npz"
    np.savez(path, **arrays)
    geo, loaded = load_geometry_archive(path)
    assert geo is None
    assert loaded._module_slices is None
    np.savez(path, **arrays, module_slices=np.array([[0, 99, 0, 6]]))
    with pytest.raises(ValueError):
        load_geometry_archive(path)
