"""ASD-POCS independently checked against dense sequential SART and TV steps."""
from unittest.mock import patch
import numpy as np
import pytest
import tigre
from tigre import algorithms as algs
from test_projection_operator import MatrixOperator, geometry


def oracle(op, data, count, beta, decay, alpha, alpha_red, rmax, epsilon, groups=None):
    x = np.zeros(4, dtype=np.float64)
    clean = np.where(op.valid_mask, data, 0).reshape(3, 2)
    residuals, steps, betas = [], [], []
    for iteration in range(count):
        old = x.copy()
        for indices in ([[0], [1], [2]] if groups is None else groups):
            a = op.matrix.reshape(3, 2, 4)[indices].reshape(-1, 4).astype(float)
            row, col = a.sum(1), a.sum(0)
            r = np.divide(clean[indices].ravel() - a@x, row, out=np.zeros_like(row), where=row>1e-6)
            update = np.divide(a.T@r, col, out=np.zeros(4), where=col>1e-6)
            x = np.maximum(x + beta*update, 0)
        dp_vec = x-old
        dp = np.linalg.norm(dp_vec)
        dd = np.linalg.norm(op.matrix@x-clean.ravel())
        if iteration == 0:
            step = alpha*dp
        old = x.copy()
        # Deterministic stand-in for TV: shrinkage exercises the adaptive
        # schedule independently of the CUDA regularizer implementation.
        x = np.maximum(x-step, 0)
        dg_vec = x-old
        dg = np.linalg.norm(dg_vec)
        steps.append(step)
        betas.append(beta)
        residuals.append(np.linalg.norm(op.matrix@x-clean.ravel()))
        if dg > rmax*dp and dd > epsilon:
            step *= alpha_red
        beta *= decay
        if (dg_vec@dp_vec / max(dg*dp, 1e-6) < -.99 and dd <= epsilon) or beta < .005:
            break
    return x, residuals, steps, betas


@pytest.mark.parametrize('name,blocksize,groups', [
    ('asd_pocs', None, [[0], [1], [2]]),
    ('awasd_pocs', None, [[0], [1], [2]]),
    ('os_awasd_pocs', 1, [[0], [1], [2]]),
    ('os_awasd_pocs', 2, [[0, 2], [1]]),
    ('os_awasd_pocs', 20, [[0, 1, 2]]),
    ('os_asd_pocs', 1, [[0], [1], [2]]),
    ('os_asd_pocs', 2, [[0, 2], [1]]),
    ('os_asd_pocs', 20, [[0, 1, 2]]),
])
def test_asd_pocs_matches_dense_oracle_masking_beta_and_tv_adaptation(name, blocksize, groups):
    op = MatrixOperator()
    data = np.arange(1, 7, dtype=np.float32).reshape(op.projection_shape)
    data[~op.valid_mask] = np.nan
    settings = dict(lmbda=.8, lmbda_red=.7, alpha=.1, alpha_red=.5, rmax=.01, maxl2err=0.)
    regularizer = '_minimize_awtv' if name in ('awasd_pocs', 'os_awasd_pocs') else '_minimize_tv'
    with patch('tigre.algorithms.operator_pocs.' + regularizer,
               side_effect=lambda x, step, *args, **kw: np.maximum(x-step, 0)) as tv:
        extra = {} if blocksize is None else dict(blocksize=blocksize)
        result, history = getattr(algs, name)(data, operator=op, niter=4, return_history=True, **extra, **settings)
    expected, residuals, steps, betas = oracle(op, data, 4, .8, .7, .1, .5, .01, 0., groups)
    np.testing.assert_allclose(result.ravel(), expected, atol=1e-6)
    for key, values in [('residual', residuals), ('tv_step', steps), ('beta', betas)]:
        np.testing.assert_allclose(history[key], values, atol=1e-6)
    assert tv.call_count == 4
    assert np.isnan(data[~op.valid_mask]).all()
    assert history['tv_step'][-1] < history['tv_step'][0]


@pytest.mark.parametrize('options', [dict(maxl2err=None), dict(maxl2err=-1),
    dict(alpha=np.nan), dict(lmbda_red=1.1), dict(alpha_red=-1), dict(tviter=1.5)])
@pytest.mark.parametrize('name', ['asd_pocs', 'os_asd_pocs', 'awasd_pocs', 'os_awasd_pocs'])
def test_reject_invalid_pocs_options(options, name):
    op = MatrixOperator()
    kwargs = dict(maxl2err=0., tviter=0)
    kwargs.update(options)
    with pytest.raises(ValueError):
        getattr(algs, name)(np.zeros(op.projection_shape), operator=op, niter=1, **kwargs)


@pytest.mark.parametrize('name', ['asd_pocs', 'os_asd_pocs', 'awasd_pocs', 'os_awasd_pocs'])
def test_zero_iterations_early_stopping_and_initializer_ownership(name):
    op = MatrixOperator()
    init = np.ones(op.volume_shape, np.float32)
    data = op.forward(init)
    result, history = getattr(algs, name)(data, operator=op, niter=0, init=init,
                                  maxl2err=0., return_history=True)
    assert not np.shares_memory(result, init) and len(history['beta']) == 0
    with patch('tigre.algorithms.operator_pocs._minimize_tv', side_effect=AssertionError('no TV step')):
        result, history = getattr(algs, name)(data, operator=op, niter=10, init=init,
            maxl2err=0., lmbda_red=0., return_history=True)
    np.testing.assert_array_equal(result, init)
    assert len(history['beta']) == 1


@pytest.mark.parametrize('blocksize', [0, -1, True, 1.5, 2.0])
@pytest.mark.parametrize('name', ['os_asd_pocs', 'os_awasd_pocs'])
def test_os_asd_rejects_invalid_blocks(blocksize, name):
    op = MatrixOperator()
    with pytest.raises(ValueError, match='blocksize'):
        getattr(algs, name)(np.zeros(op.projection_shape), operator=op,
                         niter=1, maxl2err=0., blocksize=blocksize)


@pytest.mark.parametrize('plain,ordered', [('asd_pocs', 'os_asd_pocs'), ('awasd_pocs', 'os_awasd_pocs')])
def test_single_view_blocks_match_asd_with_a_completely_missing_view(plain, ordered):
    original = MatrixOperator()
    mask = original.valid_mask.copy()
    mask[1] = False
    op = MatrixOperator(original.matrix, mask)
    data = np.ones(op.projection_shape, np.float32)
    data[~mask] = np.nan
    kwargs = dict(operator=op, niter=3, maxl2err=0., tviter=0)
    expected = getattr(algs, plain)(data, **kwargs)
    result = getattr(algs, ordered)(data, blocksize=1, **kwargs)
    np.testing.assert_array_equal(result, expected)
    assert np.isfinite(result).all()




@pytest.mark.parametrize('delta', [0., np.nan, np.inf, [1, 0, 1], [1, 2], [[1, 2, 3]], 1+1j])
@pytest.mark.parametrize('name', ['awasd_pocs', 'os_awasd_pocs'])
def test_awasd_rejects_invalid_delta_before_native_execution(delta, name):
    op = MatrixOperator()
    with pytest.raises(ValueError, match='delta'):
        getattr(algs, name)(np.zeros(op.projection_shape), operator=op, niter=0,
                        maxl2err=0., delta=delta)


@pytest.mark.parametrize('delta', [-.005])
@pytest.mark.parametrize('name', ['awasd_pocs', 'os_awasd_pocs'])
def test_awasd_native_regularization_matches_direct_kernel(delta, name):
    from _AwminTV import AwminTV
    from tigre.utilities.gpu import GpuIds
    from test_surface_projectors import fixture
    geo, surface = fixture()
    op = tigre.SurfaceOperator(geo, surface, gpuids=0)
    image = np.zeros(op.volume_shape, np.float32)
    image[2:6, 2:6, 2:6] = 1
    data = op.forward(image)
    solver = getattr(algs, name)
    options = dict(blocksize=5) if name == 'os_awasd_pocs' else {}
    data_only = solver(data, operator=op, niter=1, maxl2err=0., tviter=0, **options)
    result, history = solver(data, operator=op, niter=1, maxl2err=0.,
        tviter=2, delta=delta, return_history=True, **options)
    devices = GpuIds()
    devices.devices = [0]
    expected = AwminTV(data_only.copy(), np.float32(history['tv_step'][0]), 2, delta, devices)
    np.maximum(expected, 0, out=expected)
    np.testing.assert_allclose(result, expected, atol=2e-5, rtol=2e-5)
    wrapped = getattr(algs, name + '_surface')(data, geo, surface, niter=1,
        maxl2err=0., gpuids=0, tviter=2, delta=delta, **options)
    np.testing.assert_allclose(wrapped, result, atol=2e-5, rtol=2e-5)
    with pytest.raises(ValueError, match="directional delta"):
        solver(data, operator=op, niter=1, maxl2err=0., tviter=2, delta=[delta]*3, **options)
