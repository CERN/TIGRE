import sys
import unittest.mock

# Ensure CPU/mock compatibility if compiled CUDA extensions are not present.
for _ext in (
    "_Ax",
    "_Atb",
    "_gpuUtils",
    "_tv_proximal",
    "_minTV",
    "_minPICCS",
    "_AwminTV",
    "_RandomNumberGenerator",
):
    if _ext not in sys.modules:
        try:
            __import__(_ext)
        except ImportError:
            sys.modules[_ext] = unittest.mock.MagicMock()

import numpy as np
import pytest

from tigre.utilities.im3Dnorm import im3DNORM, l2norm


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_vector_norms_match_numpy(dtype):
    """p-norms (including 0, +-inf and NumPy integer orders) equal np.linalg.norm of the
    flattened image; the L2 norm keeps the image's dtype."""
    img = np.random.default_rng(42).standard_normal((6, 6, 6)).astype(dtype)
    flat = img.ravel()
    for order in [1, 2, 2.0, 3, -1, 0.5, 1.5, np.inf, -np.inf, np.int64(1), np.int32(2), np.float32(1.5)]:
        assert im3DNORM(img, order) == pytest.approx(float(np.linalg.norm(flat, order)), rel=1e-5)
    assert im3DNORM(img, 2).dtype == dtype
    assert im3DNORM(img, 2) == pytest.approx(float(l2norm(img)))
    assert im3DNORM(np.array([0.0, 1.5, 0.0, -2.0, 3.0]), 0) == pytest.approx(3.0)


def test_tv_norm_and_invalid_orders():
    """TV norm with forward (default), backward and central gradients; invalid input raises."""
    x = np.arange(64, dtype=np.float32).reshape(4, 4, 4)

    def tv(gz, gy, gx):
        return float(np.sum(np.sqrt(gx * gx + gy * gy + gz * gz)))

    def diffs(pad):
        return [np.pad(np.diff(x, axis=a), [pad if i == a else (0, 0) for i in range(3)]) for a in range(3)]

    forward = tv(*diffs((0, 1)))
    assert im3DNORM(x, "TV") == pytest.approx(forward)
    assert im3DNORM(x, "TV", "forward") == pytest.approx(forward)
    assert im3DNORM(x, "TV", ["forward"]) == pytest.approx(forward)
    assert im3DNORM(x, "TV", "backward") == pytest.approx(tv(*diffs((1, 0))))
    assert im3DNORM(x, "TV", "central") == pytest.approx(tv(*np.gradient(x)))

    with pytest.raises(ValueError, match="Unknown gradient option"):
        im3DNORM(x, "TV", "diagonal")
    with pytest.raises(ValueError, match="TV norm requires a 3D image"):
        im3DNORM(np.ones((4, 4)), "TV")
    for bad in [None, True, "fro", "nuc", "l2", [1, 2], np.array([1, 2])]:
        with pytest.raises(ValueError):
            im3DNORM(x, bad)
