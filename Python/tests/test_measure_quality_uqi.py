"""UQI in Measure_Quality: variance and covariance must share one normalisation.

The reference is MATLAB/Utilities/Quality_measures/UQI.m, where var() and the
covariance both divide by N-1. With a consistent estimator the index is bounded
by 1 and equals 1 exactly when the two images are identical (Wang & Bovik,
"A universal image quality index", IEEE SPL 2002).
"""

import numpy as np
import pytest

from tigre.utilities.Measure_Quality import Measure_Quality


def _uqi_matlab(real, res):
    # Line-by-line transcription of MATLAB/Utilities/Quality_measures/UQI.m
    real = real.ravel().astype(np.float64)
    res = res.ravel().astype(np.float64)
    n = real.size
    meanreal, meanres = real.mean(), res.mean()
    varreal, varres = real.var(ddof=1), res.var(ddof=1)  # MATLAB var() is N-1
    cova = np.sum((res - meanres) * (real - meanreal)) / (n - 1)
    front = (2 * cova) / (varres + varreal)
    back = (2 * meanres * meanreal) / (meanres**2 + meanreal**2)
    return front * back


@pytest.mark.parametrize("shape", [(2, 2, 2), (4, 4, 4), (16, 16, 16)])
def test_uqi_of_identical_images_is_one(shape):
    x = np.random.default_rng(0).random(shape).astype(np.float32)
    assert float(Measure_Quality(x, x, ["UQI"])) == pytest.approx(1.0, abs=1e-6)


def test_uqi_matches_matlab_reference():
    rng = np.random.default_rng(1)
    a = rng.random((8, 8, 8)).astype(np.float32)
    b = (a + 0.1 * rng.standard_normal(a.shape)).astype(np.float32)
    assert float(Measure_Quality(a, b, ["UQI"])) == pytest.approx(_uqi_matlab(a, b), rel=1e-6)


def test_uqi_never_exceeds_one():
    rng = np.random.default_rng(2)
    for _ in range(20):
        a = rng.random((3, 3, 3)).astype(np.float32)
        b = (a + 0.01 * rng.standard_normal(a.shape)).astype(np.float32)
        assert float(Measure_Quality(a, b, ["UQI"])) <= 1.0 + 1e-6
