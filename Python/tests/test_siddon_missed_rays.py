import numpy as np
import pytest
import tigre


@pytest.mark.parametrize("mode", ["cone", "parallel"])
def test_siddon_missed_rays(mode):
    """Test that Siddon projection handles missed rays properly for cone and parallel geometry."""
    geo = tigre.geometry_default(high_resolution=False)
    geo.mode = mode
    # Make detector large enough so that outer rays do not intersect the voxel volume
    geo.nDetector = np.array([256, 256])
    geo.dDetector = np.array([3.0, 3.0])
    geo.sDetector = geo.nDetector * geo.dDetector

    img = np.ones(tuple(geo.nVoxel.astype(int)), dtype=np.float32)
    angles = np.array([0.0, 0.5], dtype=np.float32)

    proj = tigre.Ax(img, geo, angles, "Siddon")

    # Missed rays must be exactly 0
    assert not np.isnan(proj).any(), f"NaNs detected in {mode} Siddon projection"
    assert (proj == 0.0).sum() > 0, f"Expected missed rays in {mode} Siddon projection"
    # Intersecting rays must be positive
    assert (proj > 0.0).sum() > 0, f"Expected intersecting rays in {mode} Siddon projection"


@pytest.mark.parametrize("mode", ["cone", "parallel"])
def test_siddon_completely_missed_volume(mode):
    """Test that Siddon projection returns all zeros when the volume is completely outside the beam."""
    geo = tigre.geometry_default(high_resolution=False)
    geo.mode = mode
    # Offset the volume far away along Y (perpendicular to optical axis)
    geo.offOrigin = np.array([0.0, 5000.0, 0.0])

    img = np.ones(tuple(geo.nVoxel.astype(int)), dtype=np.float32)
    angles = np.array([0.0], dtype=np.float32)

    proj = tigre.Ax(img, geo, angles, "Siddon")

    assert np.all(proj == 0.0), f"Expected all zeros for completely missed volume in {mode} Siddon"
