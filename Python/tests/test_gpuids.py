"""How Ax, Atb and the compiled kernels treat `gpuids`.

1. gpuids=None means "every visible GPU". Every compiled entry point converts
   it through `convert_to_c_gpuids` (_gpuUtils.pxd), which used to turn None
   into zero devices; called directly, `minTV(img, alpha, iters, None)` then
   killed the interpreter with no message. A regression here aborts the
   pytest process instead of failing an assertion.
2. Ax and Atb must not change the GpuIds they are given (#855). With fewer
   angles (Ax) or slices (Atb) than GPUs they use fewer GPUs, and they used to
   trim the caller's GpuIds to do it, so SART and friends ran on one GPU after
   their first call. Those tests replace the compiled operators by fakes and
   need no GPU.
"""
import importlib

import numpy as np
import pytest

import tigre
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb
from tigre.utilities.gpu import GpuIds

ax_mod = importlib.import_module("tigre.utilities.Ax")
atb_mod = importlib.import_module("tigre.utilities.Atb")


# ---- 1. gpuids=None (needs a GPU)
@pytest.fixture(scope="module")
def small():
    geo = tigre.geometry(mode="cone", nVoxel=np.array([32, 32, 32]), default=True)
    angles = np.linspace(0, 2 * np.pi, 30, endpoint=False, dtype=np.float32)
    geo.check_geo(angles)
    img = np.random.default_rng(0).random((32, 32, 32)).astype(np.float32)
    return geo, angles, img, Ax(img, geo, angles, "Siddon")


@pytest.mark.gpu
def test_ax_atb_accept_none(small):
    geo, angles, img, proj = small
    assert np.isfinite(Ax(img, geo, angles, "Siddon", gpuids=None)).all()
    assert np.isfinite(Atb(proj, geo, angles, backprojection_type="matched",
                           gpuids=None)).all()


@pytest.mark.gpu
def test_min_tv_accepts_none(small):
    """The case that killed the interpreter."""
    from tigre.utilities.cuda_interface._minTV import minTV
    geo, angles, img, proj = small
    out = minTV(img, 1.0, 2, None)
    assert out.shape == img.shape
    assert np.isfinite(out).all()
    assert np.any(out != 0), "denoiser returned an all-zero volume"


@pytest.mark.gpu
def test_awmin_tv_accepts_none(small):
    from tigre.utilities.cuda_interface._AwminTV import AwminTV
    geo, angles, img, proj = small
    out = AwminTV(img, 1.0, 2, -0.005, None)
    assert out.shape == img.shape
    assert np.isfinite(out).all()
    assert np.any(out != 0)


@pytest.mark.gpu
def test_tv_proximal_and_poisson_accept_none(small):
    """_tvdenoising declares the same default but is not among the built
    extensions; _tv_proximal is the one that ships."""
    from tigre.utilities.cuda_interface._tv_proximal import tvdenoise
    from tigre.utilities.cuda_interface._RandomNumberGenerator import add_poisson
    geo, angles, img, proj = small
    out = tvdenoise(img, 2, 1.0, gpuids=None)
    assert np.isfinite(out).all()
    noisy = add_poisson(np.abs(proj) + 1.0, gpuids=None)
    assert np.isfinite(noisy).all()


# ---- 2. Ax and Atb do not change the caller's GpuIds (no GPU needed)
def _two_gpus():
    return GpuIds._from_devices([0, 1])


def test_ax_trims_a_copy(monkeypatch):
    seen = []

    def fake_ax(img, geo, angles, projection_type, mode, gpuids=None):
        seen.append(list(gpuids.devices))
        return np.zeros((len(angles), *geo.nDetector), dtype=np.float32)

    monkeypatch.setattr(ax_mod, "_Ax_ext", fake_ax)
    geo = tigre.geometry(mode="cone", nVoxel=np.array([8, 8, 8]), default=True)
    img = np.zeros((8, 8, 8), dtype=np.float32)
    ids = _two_gpus()

    ax_mod.Ax(img, geo, np.array([0.0]), gpuids=ids)  # one angle: one GPU
    ax_mod.Ax(img, geo, np.linspace(0, np.pi, 4), gpuids=ids)  # four angles: both

    assert seen == [[0], [0, 1]]
    assert ids.devices == [0, 1]


def test_atb_trims_a_copy(monkeypatch):
    seen = []

    def fake_atb(proj, geo, angles, backprojection_type, mode, gpuids=None):
        seen.append(list(gpuids.devices))
        return np.zeros(geo.nVoxel, dtype=np.float32)

    monkeypatch.setattr(atb_mod, "_Atb_ext", fake_atb)
    geo = tigre.geometry(mode="cone", nVoxel=np.array([1, 8, 8]), default=True)
    angles = np.linspace(0, np.pi, 4)
    proj = np.zeros((4, *geo.nDetector), dtype=np.float32)
    ids = _two_gpus()

    atb_mod.Atb(proj, geo, angles, gpuids=ids)  # one slice: one GPU

    assert seen == [[0]]
    assert ids.devices == [0, 1]
