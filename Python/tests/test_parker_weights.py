"""Parker short-scan weights and their use in FDK filtering (CPU only, no GPU call)."""
import numpy as np
import pytest

import tigre
from tigre.utilities.filtering import filtering
from tigre.utilities.parkerweight import parkerweight


def _short_scan(nv, nu, extra=0.0, nangles=400):
    geo = tigre.geometry(mode="cone", nVoxel=np.array([8, 8, 8]), default=True)
    geo.nDetector = np.array([nv, nu])
    geo.dDetector = np.array([0.8, 0.8])
    geo.sDetector = geo.nDetector * geo.dDetector
    # fan half-angle along U
    delta = np.arctan((geo.sDetector[1] / 2 - geo.dDetector[1] / 2) / np.atleast_1d(geo.DSD)[0])
    span = np.pi + 2 * delta + extra
    angles = np.linspace(0, span + 1e-9, nangles)
    geo.check_geo(angles)
    return geo, angles


@pytest.mark.parametrize("nv, nu", [(10, 64), (64, 64), (3, 101)])
def test_weights_vary_along_u_only(nv, nu):
    geo, angles = _short_scan(nv, nu)
    proj = np.ones((len(angles), nv, nu), dtype=np.float32)
    w = parkerweight(proj, geo, angles, 1.0)
    assert w.shape == (len(angles), nv, nu)
    # every detector row gets the same weights
    np.testing.assert_array_equal(w, np.broadcast_to(w[:, :1, :], w.shape))
    assert np.all(w >= 0) and np.all(w <= 1)


@pytest.mark.parametrize("extra", [0.0, 0.3])
def test_each_ray_is_counted_over_pi(extra):
    # Parker weights make every line through the object count once: summed over
    # the scan, the weights of each detector column add up to pi radians.
    geo, angles = _short_scan(4, 128, extra=extra, nangles=2000)
    proj = np.ones((len(angles), 4, 128), dtype=np.float32)
    w = parkerweight(proj, geo, angles, 1.0)[:, 0, :]
    step = float(angles[1] - angles[0])
    np.testing.assert_allclose(w.sum(axis=0) * step, np.pi, rtol=1e-2)


def test_filtering_applies_weights_and_step_scaling():
    # With Parker weights the projections are weighted first and the angular
    # factor is 2*step per view, instead of 2*pi/nangles (as in MATLAB).
    geo, angles = _short_scan(2, 64, nangles=500)
    proj = np.random.default_rng(0).random((len(angles), 2, 64), dtype=np.float32)
    with_parker = filtering(proj.copy(), geo, angles, parker=True)
    w = parkerweight(np.ones_like(proj), geo, angles, 1.0)
    plain = filtering(proj * w, geo, angles, parker=False)
    step = float(np.mean(np.diff(angles)))
    ratio = (2 * step) / (2 * np.pi / len(angles))
    np.testing.assert_allclose(with_parker, plain * ratio, rtol=1e-5,
                               atol=1e-6 * np.abs(plain).max())


@pytest.mark.gpu
def test_fdk_short_scan_is_better_with_parker_weights():
    from tigre.algorithms import fdk
    from tigre.utilities import sample_loader

    # detector large enough that the whole head is in every projection
    geo = tigre.geometry_default(high_resolution=False)
    truth = sample_loader.load_head_phantom(geo.nVoxel).astype(np.float32)
    delta = np.arctan(geo.sDetector[1] / 2 / np.atleast_1d(geo.DSD)[0])

    def rel_error(img):
        return float(np.linalg.norm(img - truth) / np.linalg.norm(truth))

    full = np.linspace(0, 2 * np.pi, 360, endpoint=False, dtype=np.float32)
    err_full = rel_error(fdk(tigre.Ax(truth, geo, full), geo, full))

    short = np.linspace(0, np.pi + 2 * delta + 0.1, 200, dtype=np.float32)
    proj = tigre.Ax(truth, geo, short)
    err_parker = rel_error(fdk(proj.copy(), geo, short))  # on by default for short scans
    err_plain = rel_error(fdk(proj.copy(), geo, short, parker=False))
    print(f"full scan {err_full:.4f}, short scan with Parker {err_parker:.4f}, "
          f"without {err_plain:.4f}")

    assert err_parker < 0.8 * err_plain
    # a short scan with Parker weights should be close to a full scan
    assert err_parker < 2 * err_full
