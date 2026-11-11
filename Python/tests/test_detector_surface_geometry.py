"""CPU-only numerical tests of the explicit geometry module itself.

Load the file directly because importing TIGRE currently loads CUDA modules.
This is not an installed-package or native-projector test.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

path = Path(__file__).resolve().parents[1] / "tigre/utilities/detector_geometry.py"
spec = importlib.util.spec_from_file_location("surface_model_under_test", path)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
DetectorModule = model.DetectorModule
DetectorSurfaceGeometry = model.DetectorSurfaceGeometry


def test_curved_coordinates_per_view_and_missing_samples():
    source = np.array([[100.0, 0.0, 0.0], [0.0, 100.0, 0.0]])
    pixels = np.zeros((2, 2, 3, 3))
    pixels[..., 0] = [-98.0, -100.0, -98.0]
    pixels[..., 1] = [-2.0, 0.0, 2.0]
    pixels[:, 1, :, 2] = 1.0
    pixels[1, ..., 0] += 7.0
    mask = np.array([[1, 0, 1], [1, 1, 1]], dtype=bool)
    surface = DetectorSurfaceGeometry(source, pixels, mask)
    np.testing.assert_array_equal(surface.pixel_centers_mm, pixels)
    np.testing.assert_array_equal(surface.source_positions_mm, source)
    assert surface.pixel_centers_mm.dtype == np.float32
    assert surface.pixel_centers_mm.flags.c_contiguous
    assert surface.valid_mask.shape == (2, 2, 3)
    assert not surface.valid_mask[:, 0, 1].any()


def test_rigid_pose_matches_independent_quarter_turn_and_translation():
    pixels = np.array([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]])
    pose = np.array(
        [[0.0, -1.0, 0.0, 10.0], [1.0, 0.0, 0.0, 20.0], [0.0, 0.0, 1.0, 30.0], [0.0, 0.0, 0.0, 1.0]]
    )
    module = DetectorModule(pixels, pose)
    expected = np.array([[[[8.0, 21.0, 33.0], [5.0, 24.0, 36.0]]]])
    np.testing.assert_allclose(module.world_pixel_centers(2), np.repeat(expected, 2, axis=0))


def test_modules_pack_and_unpack_without_filling_gaps():
    source = [[100.0, 0.0, 0.0], [0.0, 100.0, 0.0]]
    local = np.zeros((2, 3, 3))
    first = DetectorModule(local, np.eye(4))
    poses = np.repeat(np.eye(4)[None], 2, axis=0)
    poses[1, :3, 3] = [1.0, 2.0, 3.0]
    second = DetectorModule(local, poses, column_offset=5)
    surface = DetectorSurfaceGeometry.from_modules(source, [first, second])
    assert surface.valid_mask.shape == (2, 2, 8)
    assert not surface.valid_mask[:, :, 3:5].any()
    np.testing.assert_array_equal(
        surface.pixel_centers_mm[1, :, 5:], np.broadcast_to([1.0, 2.0, 3.0], (2, 3, 3))
    )
    inputs = [np.full((2, 2, 3), 4.0), np.full((2, 2, 3), 9.0)]
    packed = surface.pack_module_data(inputs)
    assert np.all(packed[:, :, 3:5] == 0)
    for original, recovered in zip(inputs, surface.unpack_module_data(packed)):
        np.testing.assert_array_equal(original, recovered)


def test_subset_preserves_view_order_and_owns_arrays():
    surface = DetectorSurfaceGeometry(np.arange(9).reshape(3, 3), np.zeros((2, 3, 3)))
    subset = surface.subset([2, 0, 2])
    np.testing.assert_array_equal(
        subset.source_positions_mm, surface.source_positions_mm[[2, 0, 2]]
    )
    subset.pixel_centers_mm[:] = 7
    subset.valid_mask[:] = False
    assert np.all(surface.pixel_centers_mm == 0)
    assert surface.valid_mask.all()


@pytest.mark.parametrize("indices", [[], [-1], [3], [1.5], [True], [[0]]])
def test_subset_rejects_invalid_indices(indices):
    surface = DetectorSurfaceGeometry(np.zeros((3, 3)), np.zeros((2, 3, 3)))
    with pytest.raises(ValueError):
        surface.subset(indices)


@pytest.mark.parametrize(
    "source",
    [np.zeros((0, 3)), np.zeros((2, 2)), [[np.nan, 0, 0]], [[np.inf, 0, 0]], [[1e40, 0, 0]]],
)
def test_invalid_source_coordinates(source):
    with np.errstate(over="ignore"), pytest.raises(ValueError):
        DetectorSurfaceGeometry(source, np.zeros((2, 3, 3)))


@pytest.mark.parametrize("mask", [np.zeros((3, 2)), np.full((2, 3), 0.5), np.full((2, 3), np.nan)])
def test_invalid_masks(mask):
    with pytest.raises(ValueError):
        DetectorSurfaceGeometry([[100.0, 0.0, 0.0]], np.zeros((2, 3, 3)), mask)


@pytest.mark.parametrize(
    "pose", [np.diag([2.0, 1.0, 1.0, 1.0]), np.diag([-1.0, 1.0, 1.0, 1.0]), np.zeros((4, 4))]
)
def test_nonrigid_poses_are_rejected(pose):
    with pytest.raises(ValueError):
        DetectorModule(np.zeros((2, 3, 3)), pose)


def test_overlapping_layout_and_wrong_pose_count_are_rejected():
    module = DetectorModule(np.zeros((2, 3, 3)), np.eye(4))
    with pytest.raises(ValueError, match="overlap"):
        DetectorSurfaceGeometry.from_modules([[100.0, 0.0, 0.0]], [module, module])
    poses = np.repeat(np.eye(4)[None], 2, axis=0)
    with pytest.raises(ValueError, match="pose count"):
        DetectorSurfaceGeometry.from_modules(
            [[100.0, 0.0, 0.0]], [DetectorModule(np.zeros((2, 3, 3)), poses)]
        )
