"""minTV, AwminTV and minPICCS: the image boundary, and a NumPy reference.

The TV gradient uses backward differences. At the edges of the image there is
no neighbour, and the difference across the edge is 0. That was so in x and y,
but in z the GPU buffers have extra (halo) slices above and below each piece,
and at the first and last slice of the image, where those slices hold zeros,
they were read as voxels of value 0. So a constant image was not left
unchanged: its first and last slices moved. Now z is treated like x and y.

minPICCS now runs on the GPU code it shares with minTV (several GPUs, splitting,
gpuids); the reference below checks it against the formula it always computed.
"""
import numpy as np
import pytest

import tigre  # noqa: F401  (first: on Windows it puts CUDA's DLLs on the search path)

pytestmark = pytest.mark.gpu

try:
    from _minTV import minTV
    from _AwminTV import AwminTV
    from _minPICCS import minPICCS
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)


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
