"""Iterative reconstruction on curved and gapped modular detector surfaces."""
import numpy as np
import pytest

import tigre

pytest.importorskip("_AxSurface")
from tigre.algorithms.surface_algorithms import (  # noqa: E402
    cgls_surface, fast_os_sart_surface, os_sart_surface, sirt_surface)
from tigre.utilities.Ax_surface import Ax_surface  # noqa: E402
from tigre.utilities.detector_geometry import (  # noqa: E402
    DetectorModule, DetectorSurfaceGeometry)

DSO, DSD = 300.0, 500.0
NVOX = 32
# Row pitch at the isocenter should not exceed the voxel size: otherwise Siddon
# rays leave gaps between them and the problem is under-sampled (32 rows -> 0.33).
ROWS = 48
N_VIEWS = 60


def volume_geo():
    geo = tigre.geometry()
    geo.nVoxel = np.array([NVOX] * 3)
    geo.sVoxel = np.array([float(NVOX)] * 3)
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    geo.offOrigin = np.zeros(3)
    return geo


def phantom():
    ax = np.arange(NVOX) - (NVOX - 1) / 2
    z, y, x = np.meshgrid(ax, ax, ax, indexing="ij")
    img = 1.0 * ((x / 11) ** 2 + (y / 9) ** 2 + (z / 9) ** 2 < 1)
    img += 0.8 * ((x - 4) ** 2 + (y + 3) ** 2 + (z - 2) ** 2 < 9)
    img += 0.5 * ((x + 5) ** 2 + (y - 4) ** 2 + (z + 3) ** 2 < 6)
    return img.astype(np.float32)


def gantry(n_views):
    ang = np.linspace(0, 2 * np.pi, n_views, endpoint=False)
    c, s = np.cos(ang), np.sin(ang)
    rot = np.zeros((n_views, 3, 3))
    rot[:, 0, 0], rot[:, 0, 1], rot[:, 1, 0], rot[:, 1, 1], rot[:, 2, 2] = c, -s, s, c, 1
    pose = np.tile(np.eye(4), (n_views, 1, 1))
    pose[:, :3, :3] = rot
    source = rot @ np.array([DSO, 0, 0])
    return source, pose


def arc_surface():
    """Detector curved around the source, nonuniform column spacing."""
    source, pose = gantry(N_VIEWS)
    t = np.linspace(-1, 1, 64)
    gam = 0.1 * (t + 0.25 * t ** 3) / 1.25
    z = np.linspace(-30, 30, ROWS)
    local = np.zeros((ROWS, 64, 3))
    local[..., 0] = DSO - DSD * np.cos(gam)[None]
    local[..., 1] = DSD * np.sin(gam)[None]
    local[..., 2] = z[:, None]
    # DSO - DSD*cos(gam): the arc is centred on the source (at +x), curving toward -x
    module = DetectorModule(local, pose)
    return DetectorSurfaceGeometry.from_modules(source, [module])


def modular_surface():
    """Three tilted flat modules with real gaps between them."""
    source, pose = gantry(N_VIEWS)
    cols = 20
    pitch = 1.6
    modules = []
    for k, (centre, tilt) in enumerate([(-36.0, 0.14), (0.0, 0.0), (36.0, -0.14)]):
        u = pitch * (np.arange(cols) - (cols - 1) / 2)
        z = np.linspace(-30, 30, ROWS)
        local = np.zeros((ROWS, cols, 3))
        local[..., 0] = -(DSD - DSO) - u[None, :] * np.sin(tilt)
        local[..., 1] = centre + u[None, :] * np.cos(tilt)
        local[..., 2] = z[:, None]
        modules.append(DetectorModule(local, pose, column_offset=(cols + 2) * k))
    return DetectorSurfaceGeometry.from_modules(source, modules)


SURFACES = {"arc": arc_surface, "modular_gaps": modular_surface}


def rel_err(x, ref):
    return np.linalg.norm(x - ref) / np.linalg.norm(ref)


@pytest.mark.parametrize("name", list(SURFACES))
def test_cgls_reconstructs_phantom(name):
    geo = volume_geo()
    surf = SURFACES[name]()
    if name == "modular_gaps":
        assert not surf.valid_mask.all()
    truth = phantom()
    data = Ax_surface(truth, geo, surf)
    rec, res = cgls_surface(data, geo, surf, niter=100, return_residuals=True)
    assert np.all(np.diff(res) <= 1e-3 * res[0])
    # gapped modules converge slower (0.22 after 25 iterations, 0.055 after 100)
    assert rel_err(rec, truth) < 0.1, rel_err(rec, truth)


@pytest.mark.parametrize("name", list(SURFACES))
def test_sirt_reconstructs_phantom(name):
    geo = volume_geo()
    surf = SURFACES[name]()
    truth = phantom()
    data = Ax_surface(truth, geo, surf)
    rec = sirt_surface(data, geo, surf, niter=60, nonneg=True)
    assert rel_err(rec, truth) < 0.3, rel_err(rec, truth)


@pytest.mark.parametrize("name", list(SURFACES))
def test_os_sart_and_fast_variant(name):
    geo = volume_geo()
    surf = SURFACES[name]()
    truth = phantom()
    data = Ax_surface(truth, geo, surf)
    plain = os_sart_surface(data, geo, surf, niter=10, blocksize=12, nonneg=True)
    fast = fast_os_sart_surface(data, geo, surf, niter=10, blocksize=12, nonneg=True)
    # measured: arc 0.094 / 0.033, modular_gaps 0.288 / 0.217 (plain / fast)
    assert rel_err(plain, truth) < 0.32, rel_err(plain, truth)
    assert rel_err(fast, truth) < 0.9 * rel_err(plain, truth)


def test_surface_subset_matches_full_projection():
    geo = volume_geo()
    surf = modular_surface()
    truth = phantom()
    idx = np.arange(3, N_VIEWS, 7)
    full = Ax_surface(truth, geo, surf)
    part = Ax_surface(truth, geo, surf.subset(idx))
    np.testing.assert_allclose(part, full[idx], rtol=1e-5, atol=1e-4)


def test_gap_data_is_ignored():
    geo = volume_geo()
    surf = modular_surface()
    truth = phantom()
    data = Ax_surface(truth, geo, surf)
    corrupted = data.copy()
    corrupted[:, :, 20:22] = 1e3  # values in the gap columns must have no effect
    a = cgls_surface(data, geo, surf, niter=8)
    b = cgls_surface(corrupted, geo, surf, niter=8)
    np.testing.assert_allclose(a, b, rtol=1e-4, atol=1e-4)


def test_proj_shape_checked():
    geo = volume_geo()
    surf = arc_surface()
    with pytest.raises(ValueError):
        cgls_surface(np.zeros((1, 2, 3), np.float32), geo, surf)
