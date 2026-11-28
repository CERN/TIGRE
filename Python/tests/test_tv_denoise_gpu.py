"""im3ddenoise (CUDA tv_proximal) behaves like a TV proximal step (#939).

Its step sizes made the primal step negative for the first 15 iterations, so
the denoiser amplified differences (|D(x)-D(y)| up to 2.6 |x-y|) and blew up
when called repeatedly, as ossart_tv and fista do. A TV proximal step is
non-expansive and keeps the result within the range of its input.
"""
import numpy as np
import pytest

import tigre  # noqa: F401  (first: on Windows it puts CUDA's DLLs on the search path)

pytestmark = pytest.mark.gpu

try:
    from tigre.utilities.im_3d_denoise import im3ddenoise
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)


def _images():
    rng = np.random.default_rng(0)
    n = 32
    z, y, x = np.mgrid[:n, :n, :n]
    img = ((z - 16) ** 2 / 150 + (y - 16) ** 2 / 120 + (x - 16) ** 2 / 100 < 1).astype(np.float32)
    x0 = (img + 0.3 * rng.standard_normal(img.shape) - 0.2).astype(np.float32)
    x1 = (x0 + 0.01 * rng.standard_normal(x0.shape)).astype(np.float32)
    return x0, x1


@pytest.mark.parametrize("lmbda", [5, 25, 200])
def test_denoise_is_nonexpansive_and_bounded(lmbda):
    x0, x1 = _images()
    d0 = im3ddenoise(x0, 20, lmbda)
    d1 = im3ddenoise(x1, 20, lmbda)
    assert np.linalg.norm(d0 - d1) <= 1.001 * np.linalg.norm(x0 - x1)
    tol = 1e-3 * (x0.max() - x0.min())
    assert d0.min() >= x0.min() - tol and d0.max() <= x0.max() + tol


def test_repeated_denoise_does_not_grow():
    x0, _ = _images()
    z = x0.copy()
    for _ in range(30):
        z = im3ddenoise(z, 20, 25)
    assert np.isfinite(z).all()
    assert np.linalg.norm(z) <= np.linalg.norm(x0)
