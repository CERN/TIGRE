"""POCS methods whose data-consistency steps use projection operators."""

import numpy as np
from functools import partial

from .operator_algorithms import (
    _initial,
    _iterations,
    _relaxation,
    _squared_norm,
    _valid_samples,
)
from tigre.utilities.projection_operator import prepare_projections


def _finite_scalar(name, value, minimum, maximum=None):
    if (
        np.ndim(value) != 0
        or not np.isrealobj(value)
        or not np.isfinite(value)
        or value < minimum
        or (maximum is not None and value > maximum)
    ):
        raise ValueError(f'{name} must be finite in [{minimum}, {maximum or "infinity"}]')
    return float(value)


def _minimize_tv(volume, step, iterations, operator):
    """Reuse TIGRE's volume regularizer; projection backend remains separate."""
    from _minTV import minTV
    from tigre.utilities.gpu import GpuIds
    from tigre.utilities.Ax_surface import _devices

    devices = GpuIds()
    devices.devices = _devices(getattr(operator, "gpuids", None))
    return minTV(volume, np.float32(step), iterations, devices)


def _inner64(a, b):
    a, b = np.ravel(a), np.ravel(b)
    return sum(
        float(
            np.dot(a[i : i + (1 << 20)].astype(np.float64), b[i : i + (1 << 20)].astype(np.float64))
        )
        for i in range(0, a.size, 1 << 20)
    )


def _minimize_awtv(volume, step, iterations, operator, delta):
    from _AwminTV import AwminTV
    from tigre.utilities.gpu import GpuIds
    from tigre.utilities.Ax_surface import _devices

    devices = GpuIds()
    devices.devices = _devices(getattr(operator, "gpuids", None))
    if np.ndim(delta) != 0:
        raise ValueError("directional delta requires a separate native TV extension; use scalar delta")
    return AwminTV(volume, np.float32(step), iterations, delta, devices)


def awasd_pocs_operator(proj, operator, niter=20, init=None, *, delta=-0.005, **kwargs):
    """Single-view ASD-POCS with native edge-preserving weighted volume TV.

    Delta is a finite nonzero scalar or (z, y, x) vector, as in the native
    Python AwminTV binding. Both signs are accepted (the kernel squares it).
    """
    return _weighted_asd_pocs_operator(
        proj, operator, niter, init, delta=delta, blocksize=1, **kwargs
    )


def os_awasd_pocs_operator(
    proj, operator, niter=20, init=None, *, delta=-0.005, blocksize=20, **kwargs
):
    """Weighted ASD-POCS with independently normalized interleaved subsets."""
    return _weighted_asd_pocs_operator(
        proj, operator, niter, init, delta=delta, blocksize=blocksize, **kwargs
    )


def _weighted_asd_pocs_operator(proj, operator, niter, init, *, delta, blocksize, **kwargs):
    raw = np.asarray(delta)
    if raw.shape == (3,) and raw.dtype.kind in "iuf" and np.isfinite(raw).all() and np.all(raw != 0):
        raise ValueError("directional delta requires a separate native TV extension; use scalar delta")
    if raw.dtype.kind not in "iuf" or raw.shape not in ((), (3,)) or not np.isfinite(raw).all():
        raise ValueError("delta must be a finite nonzero scalar")
    values = np.array(raw, dtype=np.float32, copy=True)
    if not np.isfinite(values).all() or np.any(values == 0):
        raise ValueError("delta must be finite and nonzero in float32")
    return _asd_pocs_operator(
        proj,
        operator,
        niter,
        init,
        blocksize=blocksize,
        _regularize=partial(_minimize_awtv, delta=values),
        **kwargs,
    )


def asd_pocs_operator(proj, operator, niter=20, init=None, **kwargs):
    """Single-view ASD-POCS; see the shared kernel for options."""
    return _asd_pocs_operator(proj, operator, niter, init, blocksize=1, **kwargs)


def os_asd_pocs_operator(proj, operator, niter=20, init=None, *, blocksize=20, **kwargs):
    """ASD-POCS using interleaved subsets and subset-specific normalization."""
    return _asd_pocs_operator(proj, operator, niter, init, blocksize=blocksize, **kwargs)


def pcsd_operator(
    proj,
    operator,
    niter=20,
    init=None,
    *,
    maxl2err=None,
    lmbda=1.0,
    lmbda_red=0.99,
    tviter=20,
    nonneg=True,
    return_history=False,
):
    """PCSD with an absolute L2 tolerance and finite residual-ratio TV steps.

    Conventional PCSD is unchanged. The operator path avoids its mixed
    squared/unsquared tolerance and extra iteration. A zero first post-data
    residual disables subsequent ratio-based TV steps; an initially exact
    solution skips TV as well. It never divides by a zero reference residual.
    """
    return _asd_pocs_operator(
        proj,
        operator,
        niter,
        init,
        maxl2err=maxl2err,
        lmbda=lmbda,
        lmbda_red=lmbda_red,
        tviter=tviter,
        nonneg=nonneg,
        return_history=return_history,
        _pcsd=True,
    )


def _pcsd_keywords(kwargs):
    unknown = set(kwargs) - {"maxl2err", "lmbda", "lmbda_red", "tviter", "nonneg", "return_history"}
    if unknown:
        raise TypeError(f"Unsupported PCSD options: {sorted(unknown)}")


def os_pcsd_operator(proj, operator, niter=20, init=None, *, blocksize=20, **kwargs):
    _pcsd_keywords(kwargs)
    return _asd_pocs_operator(
        proj, operator, niter, init, blocksize=blocksize, _pcsd=True, **kwargs
    )


def aw_pcsd_operator(proj, operator, niter=20, init=None, *, delta=-0.005, **kwargs):
    _pcsd_keywords(kwargs)
    return _weighted_asd_pocs_operator(
        proj, operator, niter, init, delta=delta, blocksize=1, _pcsd=True, **kwargs
    )


def os_aw_pcsd_operator(
    proj, operator, niter=20, init=None, *, blocksize=20, delta=-0.005, **kwargs
):
    _pcsd_keywords(kwargs)
    return _weighted_asd_pocs_operator(
        proj, operator, niter, init, delta=delta, blocksize=blocksize, _pcsd=True, **kwargs
    )


def _asd_pocs_operator(
    proj,
    operator,
    niter=20,
    init=None,
    *,
    maxl2err=None,
    lmbda=1.0,
    lmbda_red=0.99,
    alpha=0.002,
    alpha_red=0.95,
    rmax=0.95,
    tviter=20,
    nonneg=True,
    return_history=False,
    blocksize=1,
    _regularize=None,
    _pcsd=False,
):
    """Shared ASD-POCS: SART/OS-SART sweeps followed by adaptive TV descent.

    ``maxl2err`` is a required absolute L2 data tolerance. There is no
    geometry-independent substitute for the legacy FDK-derived default.
    Single-view blocks use input order; larger blocks are interleaved and
    have their own masked row/column normalization. Nonnegativity is applied
    after data and TV updates.
    History records final (post-TV) residuals and the step sizes used, plus
    the pre-TV residual used in the legacy stopping/adaptation criterion.
    Conventional ``asd_pocs(data, geo, angles, niter, ...)`` is unchanged.
    """
    niter, tviter = _iterations(niter), _iterations(tviter)
    if (
        isinstance(blocksize, (bool, np.bool_))
        or not isinstance(blocksize, (int, np.integer))
        or blocksize < 1
    ):
        raise ValueError("blocksize must be a positive integer")
    if maxl2err is None:
        raise ValueError("operator ASD-POCS requires maxl2err; the legacy default depends on FDK")
    epsilon = _finite_scalar("maxl2err", maxl2err, 0)
    beta = _relaxation(lmbda)
    beta_red = _finite_scalar("lmbda_red", lmbda_red, 0, 1)
    alpha = _finite_scalar("alpha", alpha, 0)
    alpha_red = _finite_scalar("alpha_red", alpha_red, 0, 1)
    rmax = _finite_scalar("rmax", rmax, 0)
    b = prepare_projections(proj, operator)
    x = _initial(init, tuple(operator.volume_shape))
    if nonneg:
        np.maximum(x, 0, out=x)
    history = {key: [] for key in ("residual", "data_residual", "beta", "tv_step", "cosine")}
    if niter:
        ones = np.ones(operator.volume_shape, np.float32)
        subsets, weights = [], []
        nviews = operator.projection_shape[0]
        nblocks = max(1, -(-nviews // int(blocksize)))
        blocks = [np.arange(k, nviews, nblocks) for k in range(nblocks)]
        for indices in blocks:
            sub = operator.subset(indices)
            row = sub.forward(ones)
            weight = np.zeros_like(row)
            np.divide(1.0, row, out=weight, where=row > 1e-6)
            subsets.append(sub)
            weights.append(weight)
    step = 0.0
    reference_residual = None
    for iteration in range(niter):
        previous = x.copy()
        before = np.sqrt(_squared_norm(operator.forward(x) - b)) if _pcsd else None
        update_blocks = zip(blocks, subsets, weights) if not _pcsd or before > epsilon else ()
        for indices, sub, weight in update_blocks:
            # Recompute column sums to avoid caching a volume for every subset.
            col = sub.adjoint(_valid_samples(sub))
            reciprocal = np.zeros_like(col)
            np.divide(1.0, col, out=reciprocal, where=col > 1e-6)
            residual = (b[indices] - sub.forward(x)) * weight
            x += np.float32(beta) * reciprocal * sub.adjoint(residual)
            if nonneg:
                np.maximum(x, 0, out=x)
        dd = np.sqrt(_squared_norm(operator.forward(x) - b))
        dp_vector = x - previous
        dp = np.sqrt(_squared_norm(dp_vector))
        if _pcsd:
            if iteration == 0:
                step = 1.0 if before > 0 else 0.0
                reference_residual = dd
            else:
                step = before / reference_residual if reference_residual > 0 else 0.0
        elif iteration == 0:
            step = alpha * dp
        previous = x.copy()
        used_step = step
        if tviter and step > 0:
            regularize = _minimize_tv if _regularize is None else _regularize
            x = regularize(x, step, tviter, operator)
            if nonneg:
                np.maximum(x, 0, out=x)
        dg_vector = x - previous
        dg = np.sqrt(_squared_norm(dg_vector))
        if not _pcsd and dg > rmax * dp and dd > epsilon:
            step *= alpha_red
        dot = _inner64(dg_vector, dp_vector)
        cosine = dot / max(dg * dp, 1e-6)
        if return_history:
            history["residual"].append(np.sqrt(_squared_norm(operator.forward(x) - b)))
            history["data_residual"].append(dd)
            history["beta"].append(beta)
            history["tv_step"].append(used_step)
            history["cosine"].append(cosine)
        beta *= beta_red
        if (cosine < -0.99 and dd <= epsilon) or beta < 0.005:
            break
    if return_history:
        return x, {key: np.asarray(value) for key, value in history.items()}
    return x
