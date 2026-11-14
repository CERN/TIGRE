"""Thin conveniences over the same operator-aware public entry points."""

def install(namespace):
    def operator_entry(algorithm):
        def wrapped(data, operator, niter=20, **kwargs):
            return algorithm(data, operator=operator, niter=niter, **kwargs)
        return wrapped

    def surface_entry(algorithm):
        def wrapped(data, geo, surface, niter=20, gpuids=None, n_slabs=None, **kwargs):
            from tigre.utilities.projection_operator import SurfaceOperator
            return algorithm(data, operator=SurfaceOperator(geo, surface, gpuids, n_slabs),
                             niter=niter, **kwargs)
        return wrapped

    for name in ('cgls', 'sirt', 'sart', 'ossart', 'os_sart', 'fast_os_sart', 'sart_tv', 'ossart_tv', 'os_sart_tv', 'iterativereconalg', 'lsqr', 'lsmr', 'hybrid_lsqr', 'irn_tv_cgls', 'hybrid_flsqr_tv', 'ab_gmres', 'ba_gmres'):
        for suffix, factory in (("_operator", operator_entry), ("_surface", surface_entry)):
            wrapped = factory(namespace[name])
            wrapped.__name__ = name + suffix
            wrapped.__doc__ = "Explicit operator convenience for " + name
            namespace[name + suffix] = wrapped
            if name + suffix not in namespace["__all__"]:
                namespace["__all__"].append(name + suffix)
