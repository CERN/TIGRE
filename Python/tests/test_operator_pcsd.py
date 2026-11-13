"""Independent PCSD residual-gate, ratio scheduling and degenerate-case checks."""
from unittest.mock import patch
import numpy as np
import pytest
from tigre import algorithms as algs
from test_projection_operator import MatrixOperator


@pytest.mark.parametrize('epsilon', [0., 2., 100.])
def test_pcsd_dense_reference(epsilon):
    op=MatrixOperator()
    b=np.arange(1,7,dtype=np.float32).reshape(op.projection_shape)
    b[~op.valid_mask]=np.nan
    clean=np.where(op.valid_mask,b,0).ravel().astype(float)
    x=np.zeros(4); beta=.8; reference=None; residuals=[]; steps=[]; betas=[]
    for iteration in range(4):
        before=np.linalg.norm(op.matrix@x-clean)
        previous=x.copy()
        if before>epsilon:
            for view in range(3):
                a=op.matrix[2*view:2*view+2].astype(float)
                row=a.sum(1); col=a.sum(0)
                r=np.divide(clean[2*view:2*view+2]-a@x,row,out=np.zeros(2),where=row>1e-6)
                update=np.divide(a.T@r,col,out=np.zeros(4),where=col>1e-6)
                x=np.maximum(x+beta*update,0)
        dd=np.linalg.norm(op.matrix@x-clean); dpvec=x-previous
        if iteration==0:
            step=float(before>0); reference=dd
        else: step=before/reference if reference>0 else 0.
        previous=x.copy(); x=np.maximum(x-.05*step,0)
        dgvec=x-previous
        residuals.append(np.linalg.norm(op.matrix@x-clean)); steps.append(step); betas.append(beta)
        cosine=dgvec@dpvec/max(np.linalg.norm(dgvec)*np.linalg.norm(dpvec),1e-6)
        beta*=.7
        if (cosine<-.99 and dd<=epsilon) or beta<.005: break
    with patch('tigre.algorithms.operator_pocs._minimize_tv',
               side_effect=lambda x,step,*args: np.maximum(x-.05*step,0)):
        actual,h=algs.pcsd(b,operator=op,niter=4,maxl2err=epsilon,
                           lmbda=.8,lmbda_red=.7,return_history=True)
    np.testing.assert_allclose(actual.ravel(),x,atol=2e-6)
    for key,expected in [('residual',residuals),('tv_step',steps),('beta',betas)]:
        np.testing.assert_allclose(h[key],expected,atol=2e-6)
    assert np.isnan(b[~op.valid_mask]).all()


class IdentityOperator:
    volume_shape=(1,1,1)
    projection_shape=(1,1,1)
    valid_mask=None
    def forward(self,x): return x.copy()
    def adjoint(self,x): return x.copy()
    def subset(self,idx): return self


def test_absolute_tolerance_gate_and_zero_reference_are_finite():
    op=IdentityOperator(); b=np.full(op.projection_shape,.8,np.float32)
    # .8 > .7: data update must run; squaring only the residual would skip it.
    result=algs.pcsd(b,operator=op,niter=1,maxl2err=.7,tviter=0)
    np.testing.assert_array_equal(result,b)
    op.volume_shape=op.projection_shape=(1,1,2)
    b=np.array([[[.8,.4]]],np.float32)
    def perturb(x,step,*args):
        result=x.copy(); result[...,1]-=np.float32(.1*step)
        return result
    with patch('tigre.algorithms.operator_pocs._minimize_tv',
               side_effect=perturb) as tv:
        result,h=algs.pcsd(b,operator=op,niter=3,maxl2err=0.,lmbda_red=1.,return_history=True)
    assert tv.call_count==1
    np.testing.assert_array_equal(h['tv_step'],[1.,0.,0.])
    np.testing.assert_array_equal(result,b)


def test_exact_initializer_zero_iterations_and_option_validation():
    op=IdentityOperator(); init=np.ones(op.volume_shape,np.float32)
    with patch('tigre.algorithms.operator_pocs._minimize_tv',side_effect=AssertionError('unexpected TV')):
        result,h=algs.pcsd(init,operator=op,niter=3,init=init,maxl2err=0.,return_history=True)
    np.testing.assert_array_equal(result,init); assert not np.shares_memory(result,init)
    assert len(h['beta'])==3 and not np.any(h['tv_step'])
    result,h=algs.pcsd(init,operator=op,niter=0,init=init,maxl2err=0.,return_history=True)
    assert len(h['beta'])==0 and not np.shares_memory(result,init)
    for options in [{},{'maxl2err':-1},{'maxl2err':np.nan}]:
        with pytest.raises(ValueError): algs.pcsd(init,operator=op,niter=1,**options)
    with pytest.raises(TypeError): algs.pcsd(init,operator=op,niter=1,maxl2err=0.,alpha=.1)
