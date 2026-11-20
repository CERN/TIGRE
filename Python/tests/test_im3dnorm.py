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


@pytest.fixture
def sample_3d_float32():
    """Deterministic 3D float32 array for testing norms."""
    rng = np.random.default_rng(42)
    return rng.standard_normal((6, 6, 6), dtype=np.float32)


@pytest.fixture
def sample_3d_float64():
    """Deterministic 3D float64 array for testing norms."""
    rng = np.random.default_rng(42)
    return rng.standard_normal((6, 6, 6), dtype=np.float64)


def test_infinity_norms(sample_3d_float32):
    """Test positive and negative infinity norm orders."""
    img = sample_3d_float32
    flat = img.ravel()

    # np.inf
    norm_pos_inf = im3DNORM(img, np.inf)
    assert norm_pos_inf == pytest.approx(float(np.linalg.norm(flat, np.inf)))
    assert norm_pos_inf == pytest.approx(float(np.max(np.abs(flat))))

    # float('inf')
    assert im3DNORM(img, float("inf")) == pytest.approx(norm_pos_inf)

    # -np.inf
    norm_neg_inf = im3DNORM(img, -np.inf)
    assert norm_neg_inf == pytest.approx(float(np.linalg.norm(flat, -np.inf)))
    assert norm_neg_inf == pytest.approx(float(np.min(np.abs(flat))))

    # float('-inf')
    assert im3DNORM(img, float("-inf")) == pytest.approx(norm_neg_inf)


def test_integer_and_float_norm_orders(sample_3d_float64):
    """Test various valid integer and floating-point norm orders."""
    img = sample_3d_float64
    flat = img.ravel()

    for ord_val in [1, 3, 4, -1, -2, 1.5, 0.5, 2.5, np.float64(1.5), np.float32(1.5)]:
        res = im3DNORM(img, ord_val)
        expected = np.linalg.norm(flat, ord_val)
        assert res == pytest.approx(float(expected), rel=1e-6)

    # 0-norm (count of non-zero elements)
    sparse_img = np.array([0.0, 1.5, 0.0, -2.0, 3.0], dtype=np.float64)
    assert im3DNORM(sparse_img, 0) == pytest.approx(3.0)


def test_numpy_integer_types(sample_3d_float32):
    """Test NumPy integer types (np.int64, np.int32)."""
    img = sample_3d_float32
    flat = img.ravel()

    # np.int64 for non-2 norms
    res_1 = im3DNORM(img, np.int64(1))
    assert res_1 == pytest.approx(float(np.linalg.norm(flat, 1)))

    res_3 = im3DNORM(img, np.int64(3))
    assert res_3 == pytest.approx(float(np.linalg.norm(flat, 3)))

    # np.int32
    res_int32 = im3DNORM(img, np.int32(1))
    assert res_int32 == pytest.approx(float(np.linalg.norm(flat, 1)))


def test_l2_norm_consistency_and_float32_behavior(sample_3d_float32, sample_3d_float64):
    """Preserve existing optimized L2 implementation and float32 behavior."""
    img32 = sample_3d_float32

    # L2 with int 2, float 2.0, np.int64(2), np.int32(2)
    norm_int = im3DNORM(img32, 2)
    norm_float = im3DNORM(img32, 2.0)
    norm_int64 = im3DNORM(img32, np.int64(2))
    norm_int32 = im3DNORM(img32, np.int32(2))
    direct_l2 = l2norm(img32)

    # Check numerical consistency
    assert norm_int == pytest.approx(float(direct_l2))
    assert norm_float == pytest.approx(float(direct_l2))
    assert norm_int64 == pytest.approx(float(direct_l2))
    assert norm_int32 == pytest.approx(float(direct_l2))

    # Check return dtype is float32 for float32 array
    assert norm_int.dtype == np.float32
    assert norm_float.dtype == np.float32
    assert norm_int64.dtype == np.float32
    assert norm_int32.dtype == np.float32

    # Check non-float32 input (float64)
    img64 = sample_3d_float64
    norm64 = im3DNORM(img64, 2)
    assert norm64 == pytest.approx(float(np.linalg.norm(img64.ravel(), 2)))
    assert norm64.dtype == np.float64


def test_frobenius_and_nuclear_norm_raise_value_error(sample_3d_float32):
    """Frobenius and nuclear norms are invalid for vector/image raveled norms."""
    img = sample_3d_float32

    with pytest.raises(ValueError, match="Unsupported norm order"):
        im3DNORM(img, "fro")

    with pytest.raises(ValueError, match="Unsupported norm order"):
        im3DNORM(img, "nuc")


def test_unsupported_norm_values(sample_3d_float32):
    """Unsupported norm inputs must explicitly raise ValueError."""
    img = sample_3d_float32

    unsupported_inputs = [
        None,
        True,
        False,
        "unknown",
        "L2",
        "l2",
        [1, 2],
        (1, 2),
        np.array([1, 2]),
        {"norm": 2},
    ]

    for bad_order in unsupported_inputs:
        with pytest.raises(ValueError):
            im3DNORM(img, bad_order)


def test_existing_tv_norm_functionality():
    """Test TV norm calculation for 3D images with forward, backward, and central gradients."""
    x = np.arange(64, dtype=np.float32).reshape(4, 4, 4)

    # Default (forward) TV norm
    tv_default = im3DNORM(x, "TV")
    gx = np.pad(np.diff(x, axis=2), ((0, 0), (0, 0), (0, 1)), mode="constant")
    gy = np.pad(np.diff(x, axis=1), ((0, 0), (0, 1), (0, 0)), mode="constant")
    gz = np.pad(np.diff(x, axis=0), ((0, 1), (0, 0), (0, 0)), mode="constant")
    expected_forward = np.sum(np.sqrt(gx * gx + gy * gy + gz * gz))
    assert tv_default == pytest.approx(float(expected_forward))

    # Explicit 'forward'
    assert im3DNORM(x, "TV", ["forward"]) == pytest.approx(float(expected_forward))
    assert im3DNORM(x, "TV", "forward") == pytest.approx(float(expected_forward))

    # 'backward'
    tv_backward = im3DNORM(x, "TV", "backward")
    gx_b = np.pad(np.diff(x, axis=2), ((0, 0), (0, 0), (1, 0)), mode="constant")
    gy_b = np.pad(np.diff(x, axis=1), ((0, 0), (1, 0), (0, 0)), mode="constant")
    gz_b = np.pad(np.diff(x, axis=0), ((1, 0), (0, 0), (0, 0)), mode="constant")
    expected_backward = np.sum(np.sqrt(gx_b * gx_b + gy_b * gy_b + gz_b * gz_b))
    assert tv_backward == pytest.approx(float(expected_backward))

    # 'central'
    tv_central = im3DNORM(x, "TV", "central")
    gz_c, gy_c, gx_c = np.gradient(x)
    expected_central = np.sum(np.sqrt(gx_c * gx_c + gy_c * gy_c + gz_c * gz_c))
    assert tv_central == pytest.approx(float(expected_central))

    # Invalid gradient option
    with pytest.raises(ValueError, match="Unknown gradient option"):
        im3DNORM(x, "TV", "diagonal")

    # Non-3D image for TV norm
    with pytest.raises(ValueError, match="TV norm requires a 3D image"):
        im3DNORM(np.ones((4, 4)), "TV")
