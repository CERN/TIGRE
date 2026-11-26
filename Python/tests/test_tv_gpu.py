"""minTV, AwminTV and minPICCS on the GPU.

1. The image boundary. The TV gradient uses backward differences, and the
   difference across an edge of the image is 0. In z the GPU buffers have halo
   slices above and below each piece, and at the first and last slice of the
   image those were read as voxels of value 0, so a constant image was not left
   unchanged. A NumPy reference also checks minTV and minPICCS (which runs on
   minTV's GPU code since 3.3).
2. Two GPUs give the same result as one (#856). Each GPU sums its part of the
   gradient norm and the host adds the parts; that used to be in single
   precision, so the result depended on how the image was shared (up to ~3.5e-3
   of the image maximum per call). Needs 2 GPUs; skipped otherwise.
"""
import numpy as np
import pytest

import tigre  # noqa: F401  (first: on Windows it puts CUDA's DLLs on the search path)
from tigre.utilities.gpu import GpuIds

pytestmark = pytest.mark.gpu

try:
    from tigre.utilities.cuda_interface._minTV import minTV
    from tigre.utilities.cuda_interface._AwminTV import AwminTV
    from tigre.utilities.cuda_interface._minPICCS import minPICCS
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)

try:
    _N_GPUS = len(GpuIds().devices)
except Exception:  # no GPU or no driver
    _N_GPUS = 0


# ---- 1. boundary and NumPy reference
def _grad(u):
    """Backward differences (z, y, x); 0 where the neighbour is outside the image."""
    g = np.zeros((3,) + u.shape)
    g[0, 1:] = u[1:] - u[:-1]
    g[1, :, 1:] = u[:, 1:] - u[:, :-1]
    g[2, :, :, 1:] = u[:, :, 1:] - u[:, :, :-1]
    return g


def _shifted(g, axis):
    """g at the next voxel along axis (z=1, y=2, x=3 of g); 0 beyond the image."""
    out = np.zeros_like(g)
    src = [slice(None)] * 4
    dst = [slice(None)] * 4
    src[axis] = slice(1, None)
    dst[axis] = slice(None, -1)
    out[tuple(dst)] = g[tuple(src)]
    return out


def _tv_gradient(u, eps):
    g = _grad(u)
    norm = lambda h: np.sqrt((h ** 2).sum(axis=0)) + eps  # noqa: E731
    gi, gj, gk = _shifted(g, 3), _shifted(g, 2), _shifted(g, 1)
    return g.sum(axis=0) / norm(g) - gi[2] / norm(gi) - gj[1] / norm(gj) - gk[0] / norm(gk)


def _reference(img, alpha, iters, eps, prior=None, ratio=0.0):
    u = img.astype(np.float64)
    for _ in range(iters):
        d = _tv_gradient(u, eps)
        if prior is not None:
            d = (1 - ratio) * d + ratio * _tv_gradient(u - prior, eps)
        u = u - alpha * d / np.sqrt((d ** 2).sum() + 1e-7)
    return u


CALLS = {
    "minTV": lambda x: minTV(x, 0.01, 10),
    "AwminTV": lambda x: AwminTV(x, 0.01, 10, -0.005),
    "minPICCS": lambda x: minPICCS(x, np.full_like(x, 0.25), 0.01, 10, 0.5),
}


@pytest.mark.parametrize("name", sorted(CALLS))
def test_constant_image_is_left_unchanged(name):
    img = np.ones((16, 32, 24), dtype=np.float32)
    np.testing.assert_array_equal(CALLS[name](img.copy()), img)


def test_mintv_matches_reference():
    img = np.random.default_rng(0).random((12, 20, 16), dtype=np.float32)
    out = minTV(img.copy(), 0.05, 10)
    ref = _reference(img, 0.05, 10, 1e-8)
    assert np.abs(out - ref).max() < 1e-4 * np.abs(ref).max()


def test_minpiccs_matches_reference():
    rng = np.random.default_rng(1)
    img = rng.random((12, 20, 16), dtype=np.float32)
    prior = (img + 0.1 * rng.random(img.shape)).astype(np.float32)
    out = minPICCS(img.copy(), prior, 0.05, 10, 0.7)
    ref = _reference(img, 0.05, 10, 1e-6, prior=prior, ratio=0.7)
    assert np.abs(out - ref).max() < 1e-4 * np.abs(ref).max()


# ---- 2. two GPUs match one
def _gpus(devices):
    ids = GpuIds()
    ids.devices = list(devices)
    return ids


@pytest.mark.skipif(_N_GPUS < 2, reason="needs 2 GPUs")
@pytest.mark.parametrize("n", [256, 250])  # 250: the GPUs' parts do not end on reduction-block boundaries
@pytest.mark.parametrize("alpha", [0.002, 0.6])
@pytest.mark.parametrize("name", ["minTV", "AwminTV", "minPICCS"])
def test_two_gpus_match_one(name, alpha, n):
    img = np.random.default_rng(0).random((n, n, n), dtype=np.float32)
    if name == "minTV":
        call = lambda ids: minTV(img.copy(), alpha, 10, gpuids=ids)
    elif name == "AwminTV":
        call = lambda ids: AwminTV(img.copy(), alpha, 10, -0.005, gpuids=ids)
    else:  # PICCS runs on minTV's code since 3.3, so on several GPUs too
        prior = np.random.default_rng(1).random((n, n, n), dtype=np.float32)
        call = lambda ids: minPICCS(img.copy(), prior, alpha, 10, 0.5, gpuids=ids)
    one = call(_gpus([0]))
    two = call(_gpus([0, 1]))
    err = float(np.abs(two - one).max() / np.abs(one).max())
    assert err < 1e-6, err
