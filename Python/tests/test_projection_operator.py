import copy
import importlib
from unittest.mock import patch
import numpy as np
import pytest
import tigre
from tigre import algorithms as algs
from tigre.utilities.projection_operator import ProjectionOperator, SurfaceOperator, TigreOperator, prepare_projections

class MatrixOperator:
    """Positive rectangular system with an unobserved ray and voxel."""
    volume_shape = (1, 2, 2)

    def __init__(self, matrix=None, mask=None):
        self.matrix = np.array([
            [1, 2, 0, 0], [0, 1, 3, 0], [2, 0, 1, 0],
            [1, 1, 1, 0], [0, 2, 1, 0], [3, 1, 0, 0],
        ], np.float32) if matrix is None else matrix.copy()
        self.projection_shape = (len(self.matrix) // 2, 1, 2)
        self.valid_mask = np.ones(self.projection_shape, bool) if mask is None else mask.copy()
        if matrix is None:
            self.valid_mask[1, 0, 0] = False
        self.matrix[~self.valid_mask.ravel()] = 0

    def forward(self, volume):
        return (self.matrix @ volume.ravel()).reshape(self.projection_shape)

    def adjoint(self, projections):
        clean = np.where(self.valid_mask, projections, 0)
        return (self.matrix.T @ clean.ravel()).reshape(self.volume_shape)

    def subset(self, indices):
        matrix = self.matrix.reshape(self.projection_shape + (4,))[indices].reshape(-1, 4)
        return MatrixOperator(matrix, self.valid_mask[indices])


def test_cgls_matches_independent_least_squares_and_preserves_inputs():
    op = MatrixOperator()
    assert isinstance(op, ProjectionOperator)
    truth = np.array([1, 2, 3, 0], np.float32).reshape(op.volume_shape)
    b = op.forward(truth)
    b[~op.valid_mask] = np.nan
    init = np.full(op.volume_shape, .25, np.float32)
    saved = init.copy()
    result, history = algs.cgls(b, operator=op, niter=3, init=init, return_residuals=True)
    # The fourth voxel is in the nullspace: it must retain the initial value.
    expected = np.linalg.lstsq(op.matrix, np.nan_to_num(b).ravel(), rcond=None)[0]
    expected[-1] = init.ravel()[-1]
    np.testing.assert_allclose(result.ravel(), expected, atol=2e-5)
    np.testing.assert_array_equal(init, saved)
    assert np.isnan(b[~op.valid_mask]).all()
    assert history[-1] < history[0]
    assert result.dtype == np.float32


def geometry():
    geo = tigre.geometry_default(high_resolution=False)
    geo.nVoxel = np.array([3, 4, 5])
    geo.sVoxel = geo.nVoxel.astype(float)
    geo.dVoxel = np.ones(3)
    geo.nDetector = np.array([4, 5])
    geo.sDetector = geo.nDetector * geo.dDetector
    return geo


def test_conventional_nested_subset_selects_every_per_view_field_without_aliases():
    geo = geometry()
    angles = np.arange(12, dtype=float).reshape(4, 3)/10
    geo.DSO = np.arange(4, dtype=float) + 100
    geo.DSD = geo.DSO + 200
    geo.COR = np.arange(4, dtype=float)
    geo.offOrigin = np.arange(12).reshape(4, 3)
    geo.offDetector = np.arange(8).reshape(4, 2)
    geo.rotDetector = np.arange(12).reshape(4, 3)/100
    original = copy.deepcopy(geo)
    op = TigreOperator(geo, angles)
    sub = op.subset([3, 1, 3]).subset([2, 0])
    for name in ('DSO', 'DSD', 'COR', 'offOrigin', 'offDetector', 'rotDetector'):
        np.testing.assert_allclose(getattr(sub.geo, name), getattr(op.geo, name)[[3, 3]])
        np.testing.assert_array_equal(getattr(geo, name), getattr(original, name))
    np.testing.assert_array_equal(sub.angles, angles[[3, 3]])
    np.testing.assert_array_equal(sub.geo.nVoxel, [3, 4, 5])
    sub.geo.offOrigin[:] = 0
    assert op.geo.offOrigin.any()
    for invalid in [[], [True], [1.5], [-1], [4], [[0]]]:
        with pytest.raises(ValueError):
            op.subset(invalid)
    ax = importlib.import_module('tigre.utilities.Ax')
    atb = importlib.import_module('tigre.utilities.Atb')
    with patch.object(ax, 'Ax', return_value=np.zeros(op.projection_shape, np.float32)) as forward:
        op.forward(np.zeros(op.volume_shape))
        assert forward.call_args.args[0].dtype == np.float32
        assert forward.call_args.args[3] == 'Siddon'
    with patch.object(atb, 'Atb', return_value=np.zeros(op.volume_shape, np.float32)) as adjoint:
        op.adjoint(np.zeros(op.projection_shape))
        assert adjoint.call_args.args[3] == 'matched'


@pytest.mark.parametrize("solver", [algs.cgls_operator, algs.sirt_operator])
def test_initialization_validation_and_exact_solution(solver):
    op = MatrixOperator()
    init = np.ones(op.volume_shape, np.float32)
    data = op.forward(init)
    np.testing.assert_array_equal(solver(data, op, niter=0, init=init), init)
    np.testing.assert_allclose(solver(data, op, niter=2, init=init), init)
    for count in [-1, True, 1.5]:
        with pytest.raises(ValueError, match="niter"):
            solver(data, op, niter=count)
    with pytest.raises(ValueError):
        solver(data, op, init="FDK")
    with pytest.raises(ValueError):
        solver(data, op, init=np.full(op.volume_shape, np.nan))

def test_sirt_dense_update():
    op = MatrixOperator()
    data = np.arange(1, 7, dtype=np.float32).reshape(op.projection_shape)
    data[~op.valid_mask] = np.nan
    a = op.matrix.astype(float)
    rows, cols = a.sum(1), a.sum(0)
    residual = np.divide(np.nan_to_num(data).ravel(), rows, out=np.zeros(6), where=rows>0)
    expected = np.divide(a.T @ residual, cols, out=np.zeros(4), where=cols>0)
    np.testing.assert_allclose(algs.sirt(data, operator=op, niter=1).ravel(), expected, atol=1e-6)

def test_legacy_dispatch_is_preserved():
    op = MatrixOperator()
    with pytest.raises(ValueError, match="omit geo"):
        algs.cgls(np.zeros(op.projection_shape), object(), None, 1, operator=op)
    module = importlib.import_module("tigre.algorithms.art_family_algorithms")
    with patch.object(module.SIRT, "__init__", return_value=None) as init, \
         patch.object(module.SIRT, "run_main_iter"), \
         patch.object(module.SIRT, "getres", return_value="legacy"):
        module.SIRT.computel2 = False
        module.SIRT.Quameasopts = None
        geo = type("Geometry", (), {"check_geo": lambda self, angles: None})()
        assert algs.sirt(None, geo, np.zeros(2), 3, verbose=False) == "legacy"
        assert init.call_args.kwargs == {"verbose": False}
        del module.SIRT.computel2, module.SIRT.Quameasopts
