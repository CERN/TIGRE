"""Shared CGLS/SIRT kernels for conventional and surface adapters."""
import numpy as np
from tigre.utilities.projection_operator import prepare_projections, _array

def reconstruct(name, proj, operator, niter, _algorithm_class=None, **kwargs):
    name = name.lower() if isinstance(name, str) else "iterativereconalg"
    name = {"ossart": "os_sart", "ossart_tv": "os_sart_tv"}.get(name, name)
    if name in ("asd_pocs", "os_asd_pocs", "awasd_pocs", "os_awasd_pocs", "pcsd", "os_pcsd", "aw_pcsd", "os_aw_pcsd"):
        from . import operator_pocs
        return getattr(operator_pocs, name + "_operator")(proj, operator, niter=niter, **kwargs)
    solvers = {'cgls': cgls_operator, 'sirt': sirt_operator, 'os_sart': os_sart_operator, 'fast_os_sart': fast_os_sart_operator}
    if name not in solvers:
        if name not in ('sart', 'sart_tv', 'os_sart_tv', 'iterativereconalg', 'lsqr', 'lsmr', 'hybrid_lsqr', 'irn_tv_cgls', 'hybrid_flsqr_tv', 'ab_gmres', 'ba_gmres', 'fista', 'ista', 'mlem'):
            raise NotImplementedError(f"{name} is not included in this operator slice")
        from .operator_recon_adapter import run_existing
        return run_existing(_algorithm_class, proj, operator, niter, **kwargs)
    if "noneg" in kwargs:
        if "nonneg" in kwargs:
            raise TypeError("specify only one of noneg and nonneg")
        kwargs["nonneg"] = kwargs.pop("noneg")
    return solvers[name](proj, operator, niter=niter, **kwargs)

def _iterations(niter):
    if isinstance(niter, (bool, np.bool_)) or not isinstance(niter, (int, np.integer)) or niter < 0:
        raise ValueError("niter must be a non-negative integer")
    return int(niter)


def _initial(init, shape):
    if isinstance(init, str):
        raise ValueError(
            "operator solvers accept only a volume initializer or None; FDK is geometry-specific"
        )
    return np.zeros(shape, np.float32) if init is None else _array(init, shape, "init")


def _valid_samples(operator):
    return (
        np.ones(operator.projection_shape, np.float32)
        if operator.valid_mask is None
        else np.asarray(operator.valid_mask, dtype=np.float32)
    )


def _relaxation(value):
    if np.ndim(value) != 0 or not np.isfinite(value) or value <= 0:
        raise ValueError("lmbda must be finite and positive")
    return value


def _squared_norm(value):
    """Float64 accumulation with bounded scratch space for large CT arrays."""
    flat = np.ravel(value)
    total = 0.0
    for start in range(0, flat.size, 1 << 20):
        chunk = flat[start : start + (1 << 20)].astype(np.float64)
        total += float(np.dot(chunk, chunk))
    return total


def cgls_operator(proj, operator, niter=20, init=None, return_residuals=False):
    """Conjugate gradient least squares with a projection operator.

    Backend settings belong to ``operator``.
    """
    niter = _iterations(niter)
    b = prepare_projections(proj, operator)
    shape = tuple(operator.volume_shape)
    x = _initial(init, shape)
    if niter == 0:
        return (x, np.empty(0)) if return_residuals else x
    r = b - (operator.forward(x) if init is not None else 0)
    p = operator.adjoint(r.astype(np.float32))
    gamma = _squared_norm(p)
    residuals = []
    for _ in range(int(niter)):
        if gamma == 0.0:
            break
        q = operator.forward(p)
        qq = _squared_norm(q)
        if qq == 0.0:
            break
        alpha = gamma / qq
        x += np.float32(alpha) * p
        r = r - np.float32(alpha) * q
        residuals.append(float(np.sqrt(_squared_norm(r))))
        s = operator.adjoint(r.astype(np.float32))
        gamma_new = _squared_norm(s)
        p = s + np.float32(gamma_new / gamma) * p
        gamma = gamma_new
    return (x, np.array(residuals)) if return_residuals else x


def sirt_operator(proj, operator, niter=50, init=None, lmbda=1.0, nonneg=False):
    """SIRT with row/column sums computed using the selected operator."""
    niter = _iterations(niter)
    b = prepare_projections(proj, operator)
    lmbda = _relaxation(lmbda)
    shape = tuple(operator.volume_shape)
    x = _initial(init, shape)
    if niter == 0:
        return x

    ones_vol = np.ones(shape, np.float32)
    row_sum = operator.forward(ones_vol)
    w = np.zeros_like(row_sum)
    np.divide(1.0, row_sum, out=w, where=row_sum > 1e-6)
    col_sum = operator.adjoint(_valid_samples(operator))
    v = np.zeros_like(col_sum)
    np.divide(1.0, col_sum, out=v, where=col_sum > 1e-6)

    for _ in range(int(niter)):
        resid = (b - operator.forward(x)) * w
        x += np.float32(lmbda) * v * operator.adjoint(resid.astype(np.float32))
        if nonneg:
            np.maximum(x, 0, out=x)
    return x

def os_sart_operator(
    proj, operator, niter=10, blocksize=20, init=None, lmbda=1.0, nonneg=False, momentum=False
):
    """Ordered-subset SART with a projection operator.

    Views are split into interleaved blocks. With ``momentum=True`` the iterate
    is extrapolated before every pass with the Nesterov schedule used by
    ``tigre.algorithms.fast_os_sart`` (an acceleration heuristic, not a
    convergent method). The per-block back-projection normalisation is
    recomputed each time to keep memory at one extra volume. Backend settings belong to ``operator``.
    """
    niter = _iterations(niter)
    b = prepare_projections(proj, operator)
    if (
        isinstance(blocksize, (bool, np.bool_))
        or not isinstance(blocksize, (int, np.integer))
        or blocksize < 1
    ):
        raise ValueError(f"blocksize must be a positive integer, got {blocksize!r}")
    lmbda = _relaxation(lmbda)
    shape = tuple(operator.volume_shape)
    x = _initial(init, shape)
    if niter == 0:
        return x
    n_views = b.shape[0]
    n_blocks = max(1, -(-n_views // int(blocksize)))
    blocks = [np.arange(k, n_views, n_blocks) for k in range(n_blocks)]
    ones_vol = np.ones(shape, np.float32)

    subs, weights = [], []
    for idx in blocks:
        sub = operator.subset(idx)
        row_sum = sub.forward(ones_vol)
        w = np.zeros_like(row_sum)
        np.divide(1.0, row_sum, out=w, where=row_sum > 1e-6)
        subs.append(sub)
        weights.append(w)

    t = 1.0
    x_old = x.copy()
    for _ in range(int(niter)):
        if momentum:
            t_old, t = t, (1.0 + np.sqrt(1.0 + 4.0 * t * t)) / 2.0
            beta = np.float32((t_old - 1.0) / t)
            for k in range(x.shape[0]):
                saved = x[k].copy()
                x[k] += beta * (x[k] - x_old[k])
                x_old[k] = saved
        for idx, sub, w in zip(blocks, subs, weights):
            resid = (b[idx] - sub.forward(x)) * w
            col_sum = sub.adjoint(_valid_samples(sub))
            v = np.zeros_like(col_sum)
            np.divide(1.0, col_sum, out=v, where=col_sum > 1e-6)
            x += np.float32(lmbda) * v * sub.adjoint(resid.astype(np.float32))
            if nonneg:
                np.maximum(x, 0, out=x)
    return x


def fast_os_sart_operator(proj, operator, niter=10, **kwargs):
    """OS-SART with Nesterov momentum with a projection operator."""
    return os_sart_operator(proj, operator, niter=niter, momentum=True, **kwargs)
