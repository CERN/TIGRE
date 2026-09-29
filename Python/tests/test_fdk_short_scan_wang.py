"""FDK must not apply Wang detector-offset weights on a short scan.

Wang weighting (displaced-detector redundancy weights) assumes a full circle:
it ramps one side of the detector down and relies on the opposing (beta + pi)
views to restore uniform coverage. On a short scan those views do not exist,
so the ramp survives as a one-sided "teardrop" shading. Any non-zero
offDetector - a sub-pixel calibration value is enough - switches the weights
on, so short scans with a calibrated detector offset were silently wrong.

The classification tests are pure NumPy. The reconstruction test needs a GPU
and is skipped without one.
"""
import numpy as np
import pytest

from tigre.algorithms.single_pass_algorithms import is_short_scan


class TestIsShortScan:
    def test_full_circle_endpoint_false_is_not_short(self):
        assert not is_short_scan(np.linspace(0, 2 * np.pi, 360, endpoint=False))

    def test_full_circle_endpoint_true_is_not_short(self):
        assert not is_short_scan(np.linspace(0, 2 * np.pi, 361))

    def test_200_degrees_is_short(self):
        assert is_short_scan(np.linspace(0, np.deg2rad(200), 200))

    def test_half_circle_is_short(self):
        assert is_short_scan(np.linspace(0, np.pi, 180, endpoint=False))

    def test_euler_triplets_use_first_column(self):
        a = np.linspace(0, 2 * np.pi, 360, endpoint=False)
        full = np.stack([a, np.zeros_like(a), np.zeros_like(a)], axis=1)
        short = np.stack([a[:200], np.zeros(200), np.zeros(200)], axis=1)
        assert not is_short_scan(full)
        assert is_short_scan(short)

    def test_descending_angles(self):
        assert not is_short_scan(np.linspace(2 * np.pi, 0, 360, endpoint=False))
        assert is_short_scan(np.linspace(np.deg2rad(200), 0, 200))

    def test_single_angle_is_short(self):
        assert is_short_scan(np.array([0.3]))


def _gpu_available():
    try:
        from tigre.utilities.gpu import getGpuNames
        return len(getGpuNames()) > 0
    except Exception:
        return False


@pytest.mark.skipif(not _gpu_available(), reason="needs a CUDA GPU")
def test_short_scan_with_subpixel_offset_has_no_teardrop():
    """A 200-degree scan of a uniform cylinder with a 0.1 mm (sub-pixel)
    detector offset: the interior must reconstruct left/right balanced, as
    the centred short scan does. Before the gate the Wang ramp made one half
    much brighter than the other."""
    import tigre
    import tigre.algorithms as algs

    geo = tigre.geometry_default(high_resolution=False)
    geo.nVoxel = np.array([64, 64, 64])
    geo.sVoxel = np.array([64.0, 64.0, 64.0])
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    zz, yy, xx = np.mgrid[:64, :64, :64]
    phantom = (((yy - 31.5) ** 2 + (xx - 31.5) ** 2 <= 20 ** 2) * 0.02).astype(np.float32)
    core = (yy - 31.5) ** 2 + (xx - 31.5) ** 2 <= 14 ** 2
    angles = np.linspace(0, np.deg2rad(200), 200)
    proj = tigre.Ax(phantom, geo, angles)

    def lr_ratio(vol):
        left = float(vol[16:48, :, :32][core[16:48, :, :32]].mean())
        right = float(vol[16:48, :, 32:][core[16:48, :, 32:]].mean())
        return left / right

    ref = lr_ratio(algs.fdk(proj, geo, angles))                  # centred
    geo.offDetector = np.array([0.0, 0.1])                       # sub-pixel
    gated = lr_ratio(algs.fdk(proj, geo, angles))
    assert gated == pytest.approx(ref, rel=0.05)
