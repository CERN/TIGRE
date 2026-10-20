"""Tests for the threaded FDK/FBP filtering step (CPU only, no GPU call)."""
import numpy as np
import pytest

import tigre
import tigre.utilities.filtering as filtering_module
from tigre.utilities.filtering import default_workers, filter as make_filter, filtering, nextpow2
from tigre.utilities.filtering import ramp_flat


def _geo(nangles, nv, nu, filt=None):
    geo = tigre.geometry(mode="cone", nVoxel=np.array([8, 8, 8]), default=True)
    geo.nDetector = np.array([nv, nu])
    geo.dDetector = np.array([0.8, 0.8])
    geo.sDetector = geo.nDetector * geo.dDetector
    angles = np.linspace(0, 2 * np.pi, nangles, endpoint=False, dtype=np.float32)
    geo.check_geo(angles)
    geo.filter = filt
    return geo, angles


def _proj(nangles, nv, nu):
    return np.random.default_rng(nangles + nv + nu).random((nangles, nv, nu), dtype=np.float32)


def _plain_filter(proj, geo, angles):
    """One projection at a time, float64, no packing, no blocks, no threads."""
    nu = int(geo.nDetector[1])
    filt_len = max(64, 2 ** nextpow2(2 * nu))
    filt = make_filter(geo.filter, ramp_flat(filt_len)[0], filt_len, 1).astype(np.float64)
    padding = (filt_len - nu) // 2
    scale = (geo.DSD[0] / geo.DSO[0]) * (2 * np.pi / len(angles)) / (4 * geo.dDetector[1])
    out = np.empty(proj.shape, dtype=np.float64)
    for i in range(proj.shape[0]):
        padded = np.zeros((proj.shape[1], filt_len))
        padded[:, padding:padding + nu] = proj[i]
        filtered = np.fft.ifft(np.fft.fft(padded, axis=1) * filt, axis=1).real
        out[i] = filtered[:, padding:padding + nu] * scale
    return out


@pytest.mark.parametrize("shape", [(6, 40, 64), (7, 40, 64), (1, 3, 50), (9, 1, 200), (4, 300, 256)])
@pytest.mark.parametrize("filt", [None, "ram_lak", "shepp_logan", "cosine", "hamming", "hann"])
def test_matches_a_plain_per_projection_filter(shape, filt):
    geo, angles = _geo(*shape, filt=filt)
    proj = _proj(*shape)
    expected = _plain_filter(proj, geo, angles)
    out = filtering(proj.copy(), geo, angles, parker=False)
    assert out.dtype == np.float32
    np.testing.assert_allclose(out, expected, rtol=0, atol=2e-6 * np.abs(expected).max())


@pytest.mark.parametrize("shape", [(6, 40, 64), (7, 13, 33), (9, 1, 200), (5, 300, 256)])
def test_same_bits_whatever_the_number_of_workers(shape, monkeypatch):
    geo, angles = _geo(*shape)
    proj = _proj(*shape)
    serial = filtering(proj.copy(), geo, angles, parker=False, workers=1)
    for workers in (None, 2, 3, 8):
        assert np.array_equal(serial, filtering(proj.copy(), geo, angles, parker=False,
                                                workers=workers))
    # the smallest blocks the function will use: many blocks per projection
    monkeypatch.setattr(filtering_module, "BLOCK_BYTES", 1)
    for workers in (1, 4):
        assert np.array_equal(serial, filtering(proj.copy(), geo, angles, parker=False,
                                                workers=workers))


def test_filters_in_place_and_returns_its_argument():
    geo, angles = _geo(4, 20, 32)
    proj = _proj(4, 20, 32)
    before = proj.copy()
    out = filtering(proj, geo, angles, parker=False)
    assert out is proj
    assert not np.array_equal(before, proj)


def test_error_in_a_worker_thread_reaches_the_caller():
    geo, angles = _geo(8, 600, 512)
    with pytest.raises(ValueError):
        filtering(_proj(8, 600, 500), geo, angles, parker=False, workers=4)


def test_default_workers(monkeypatch):
    assert 1 <= default_workers() <= filtering_module.MAX_WORKERS
    monkeypatch.setattr(filtering_module, "MAX_WORKERS", 1)
    assert default_workers() == 1
