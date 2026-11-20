"""CTnoise.add: a new noise realisation every call, unless a seed is given.

The GPU random number generator used to be seeded with 1234 on every call, so
every call added exactly the same noise.
"""
import numpy as np
import pytest

import tigre  # noqa: F401  (first: on Windows it puts CUDA's DLLs on the search path)

pytestmark = pytest.mark.gpu

try:
    from tigre.utilities import CTnoise
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)


def _proj():
    return np.random.default_rng(0).random((4, 16, 24), dtype=np.float32) * 5


def test_two_calls_give_different_noise():
    p = _proj()
    assert not np.array_equal(CTnoise.add(p), CTnoise.add(p))


def test_same_seed_same_noise():
    p = _proj()
    a = CTnoise.add(p, seed=7)
    np.testing.assert_array_equal(CTnoise.add(p, seed=7), a)
    assert not np.array_equal(CTnoise.add(p, seed=8), a)


def test_numpy_seed_makes_it_reproducible():
    p = _proj()
    np.random.seed(3)
    a = CTnoise.add(p)
    np.random.seed(3)
    np.testing.assert_array_equal(CTnoise.add(p), a)


def test_negative_seed_is_refused():
    with pytest.raises(ValueError):
        CTnoise.add(_proj(), seed=-1)
