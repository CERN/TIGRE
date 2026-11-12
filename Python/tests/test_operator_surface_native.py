import numpy as np
import pytest
from tigre import algorithms as algs
from tigre.utilities.projection_operator import SurfaceOperator
from tigre.utilities.gpu import getGpuNames
from test_surface_projectors import fixture

pytestmark = pytest.mark.skipif(not getGpuNames(), reason="CUDA GPU required")

def test_surface_adapter_subsets_adjoint_and_solver_wrappers():
    geo, surface = fixture()
    op = SurfaceOperator(geo, surface, gpuids=0, n_slabs=2)
    rng = np.random.default_rng(7)
    x = rng.random(op.volume_shape, dtype=np.float32)
    data = rng.random(op.projection_shape, dtype=np.float32)
    data[~op.valid_mask] = np.nan
    clean = np.where(op.valid_mask, data, 0)
    np.testing.assert_allclose(np.vdot(op.forward(x), clean), np.vdot(x, op.adjoint(data)), rtol=2e-5)
    np.testing.assert_allclose(op.subset([2, 0, 2]).subset([1]).forward(x), op.forward(x)[[0]], atol=1e-5)
    projections = op.forward(x)
    for name in ("cgls", "sirt"):
        result = getattr(algs, name)(projections, operator=op, niter=2)
        wrapped = getattr(algs, name + "_surface")(projections, geo, surface, niter=2, gpuids=0, n_slabs=2)
        np.testing.assert_allclose(result, wrapped, atol=2e-5)
        assert np.linalg.norm(op.forward(result)-projections) < np.linalg.norm(projections)
