import numpy as np
import pytest
from tigre import algorithms as algs
from test_projection_operator import MatrixOperator

@pytest.mark.parametrize('name',['fista','ista'])
def test_proximal_data_steps_against_dense_gradient(name):
    op=MatrixOperator(); b=op.forward(np.ones(op.volume_shape,np.float32)); b[~op.valid_mask]=np.nan
    L=200.; x=np.zeros(4); previous=x.copy(); t=1.
    clean=np.nan_to_num(b).ravel()
    for _ in range(2):
        next=np.maximum(x+2/L*op.matrix.T@(clean-op.matrix@x),0)
        if name=='fista':
            new_t=(1+np.sqrt(1+4*t*t))/2
            x=next+(t-1)/new_t*(next-previous); previous=next; t=new_t
        else: x=next
    result=getattr(algs,name)(b,operator=op,niter=2,hyper=L,tvlambda=0.,tviter=0,verbose=False)
    np.testing.assert_allclose(result.ravel(),x,atol=1e-6)


def test_mlem_multiplicative_update_and_masked_sensitivity():
    op=MatrixOperator(); b=op.forward(np.arange(1,5,dtype=np.float32).reshape(op.volume_shape))
    b[~op.valid_mask]=np.nan
    den=op.matrix@np.ones(4); ratio=np.divide(np.nan_to_num(b).ravel(),den,out=np.zeros(6),where=den>1e-6)
    sensitivity=op.matrix.sum(0)
    expected=np.divide(op.matrix.T@np.clip(ratio,0,10),sensitivity,out=np.zeros(4),where=sensitivity>0)
    result=algs.mlem(b,operator=op,niter=1,verbose=False)
    np.testing.assert_allclose(result.ravel(),expected,atol=1e-6)


def test_zero_iterations_and_invalid_algorithm_specific_settings():
    op=MatrixOperator(); data=np.zeros(op.projection_shape,np.float32)
    np.testing.assert_array_equal(algs.mlem(data,operator=op,niter=0),1)
    for name in ('lsqr','lsmr','sart','fista','mlem'):
        result,errors=getattr(algs,name)(data,operator=op,niter=0,computel2=True)
        assert result.shape==op.volume_shape and errors.shape==(1,0)
    with pytest.raises(ValueError): algs.os_sart(np.ones_like(data),operator=op,niter=1,blocksize=0)
    with pytest.raises(ValueError): algs.fista(np.ones_like(data),operator=op,niter=1,hyper=-1)
    with pytest.raises(ValueError): algs.mlem(-np.ones_like(data),operator=op,niter=1)
