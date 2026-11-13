import numpy as np
import pytest
from tigre import algorithms as algs
from test_projection_operator import MatrixOperator

@pytest.mark.parametrize('name',['lsqr','lsmr','ab_gmres','ba_gmres'])
def test_krylov_solution_against_dense_least_squares(name):
    op=MatrixOperator(); truth=np.array([1,2,3,0],np.float32).reshape(op.volume_shape)
    b=op.forward(truth); b[~op.valid_mask]=np.nan
    options=dict(operator=op,niter=3,verbose=False)
    if name=='lsmr': options['lmbda']=0.
    actual=getattr(algs,name)(b,**options)
    expected=np.linalg.lstsq(op.matrix,np.nan_to_num(b).ravel(),rcond=None)[0]
    np.testing.assert_allclose(actual.ravel(),expected,atol=1e-4,rtol=1e-4)


def test_lsmr_damping_against_augmented_dense_system():
    op=MatrixOperator(); b=op.forward(np.ones(op.volume_shape,np.float32))
    damp=.7
    result=algs.lsmr(b,operator=op,niter=3,lmbda=damp,verbose=False)
    augmented=np.concatenate([op.matrix,damp*np.eye(4)])
    expected=np.linalg.lstsq(augmented,np.r_[b.ravel(),np.zeros(4)],rcond=None)[0]
    np.testing.assert_allclose(result.ravel(),expected,atol=1e-4,rtol=1e-4)


@pytest.mark.parametrize("name", ["lsqr", "lsmr", "hybrid_lsqr", "ab_gmres", "ba_gmres", "irn_tv_cgls", "hybrid_flsqr_tv"])
def test_exact_initializers_and_invalid_options(name):
    op=MatrixOperator()
    x=np.ones(op.volume_shape,np.float32)
    options={"niter_outer":1} if name=="irn_tv_cgls" else {}
    actual=getattr(algs,name)(op.forward(x),operator=op,niter=0,init=x,**options)
    np.testing.assert_array_equal(actual,x)
    with pytest.raises(TypeError):
        getattr(algs,name)(op.forward(x),operator=op,niter=0,UnknownOption=1)
