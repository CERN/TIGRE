"""Shared CGLS/SIRT kernels for conventional and surface adapters."""
import numpy as np
from tigre.utilities.projection_operator import prepare_projections, _array

def reconstruct(name, proj, operator, niter, _algorithm_class=None, **kwargs):
    name = name.lower() if isinstance(name, str) else "iterativereconalg"
    name = {"ossart": "os_sart", "ossart_tv": "os_sart_tv"}.get(name, name)
    solvers = {'cgls': cgls_operator, 'sirt': sirt_operator}
    if name not in solvers:
        raise NotImplementedError(f"{name} does not yet support projection operators")
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
