import numpy as np
import pytest
from tigre import algorithms as algs
from test_projection_operator import MatrixOperator

@pytest.mark.parametrize('solver', [algs.sirt_operator, algs.os_sart_operator, algs.fast_os_sart_operator])
def test_normalization_and_interleaved_subsets_match_matrix_oracle(solver):
    op = MatrixOperator()
    b = np.arange(1, 7, dtype=np.float32).reshape(op.projection_shape)
    b[~op.valid_mask] = np.inf
    is_sirt = solver is algs.sirt_operator
    options = {} if is_sirt else dict(blocksize=2)
    actual = solver(b, op, niter=3, lmbda=.8, nonneg=True, **options)
    expected = np.zeros(4)
    previous = expected.copy()
    t = 1.
    groups = [np.arange(3)] if is_sirt else [np.array([0, 2]), np.array([1])]
    for _ in range(3):
        if solver is algs.fast_os_sart_operator:
            new_t = (1 + np.sqrt(1 + 4*t*t))/2
            saved = expected.copy()
            expected += (t-1)/new_t * (expected-previous)
            previous, t = saved, new_t
        for indices in groups:
            a = op.matrix.reshape(3, 2, 4)[indices].reshape(-1, 4).astype(float)
            data = np.where(op.valid_mask[indices], b[indices], 0).ravel()
            row, col = a.sum(axis=1), a.sum(axis=0)
            residual = np.divide(data - a@expected, row, out=np.zeros_like(row), where=row>1e-6)
            correction = np.divide(a.T@residual, col, out=np.zeros_like(col), where=col>1e-6)
            expected = np.maximum(expected + .8*correction, 0)
    np.testing.assert_allclose(actual.ravel(), expected, atol=2e-6)


@pytest.mark.parametrize("name", ["sart", "ossart", "os_sart", "sart_tv", "ossart_tv", "os_sart_tv", "iterativereconalg"])
def test_algebraic_adaptation(name):
    op=MatrixOperator()
    data=op.forward(np.ones(op.volume_shape,np.float32))
    options={"tviter":0} if "tv" in name else {}
    actual=getattr(algs,name)(data,operator=op,niter=2,**options)
    assert np.isfinite(actual).all()
    assert np.linalg.norm(op.forward(actual)-data)<np.linalg.norm(data)
    np.testing.assert_allclose(actual,getattr(algs,name+"_operator")(data,op,niter=2,**options))
