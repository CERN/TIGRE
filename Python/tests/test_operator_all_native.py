"""Native smoke and wrapper equivalence against the submitted projector API."""
import numpy as np
import pytest
from tigre import algorithms as algs
from tigre.utilities.projection_operator import SurfaceOperator
from test_surface_projectors import fixture

@pytest.mark.parametrize("name", ['cgls', 'sirt', 'sart', 'os_sart', 'fast_os_sart', 'sart_tv', 'os_sart_tv', 'iterativereconalg', 'lsqr', 'lsmr', 'hybrid_lsqr', 'irn_tv_cgls', 'hybrid_flsqr_tv', 'ab_gmres', 'ba_gmres', 'fista', 'ista', 'mlem', 'asd_pocs', 'os_asd_pocs', 'awasd_pocs', 'os_awasd_pocs', 'pcsd', 'os_pcsd', 'aw_pcsd', 'os_aw_pcsd'])
def test_native_algorithm_and_surface_wrapper(name):
    geo,surface=fixture()
    op=SurfaceOperator(geo,surface,gpuids=0)
    x=np.zeros(op.volume_shape,np.float32)
    x[1:3,1:4,1:5]=1
    x[2,2,3]=1.7
    data=op.forward(x)
    data[~op.valid_mask]=np.nan
    options={}
    if name in ("sart_tv","os_sart_tv"):
        options.update(tviter=2,tvlambda=50.)
    if name in ("fista","ista"):
        options.update(tviter=2,tvlambda=.1)
    if "pocs" in name or "pcsd" in name:
        options.update(maxl2err=0.,tviter=2)
    if name=="irn_tv_cgls":
        options["niter_outer"]=1
    result=getattr(algs,name)(data,operator=op,niter=2,**options)
    wrapped=getattr(algs,name+"_surface")(data,geo,surface,niter=2,gpuids=0,**options)
    assert result.dtype==np.float32 and np.isfinite(result).all()
    np.testing.assert_allclose(result,wrapped,rtol=2e-3,atol=2e-4)
    residual=op.forward(result)-np.where(op.valid_mask,data,0)
    assert np.linalg.norm(residual)<1.1*np.linalg.norm(np.nan_to_num(data))
