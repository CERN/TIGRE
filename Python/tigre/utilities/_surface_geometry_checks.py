"""CPU geometry diagnostics only; no reconstruction or accuracy certification."""
import numpy as np

class _Fan:
    """Per-column / per-row tables of a validated circular-orbit surface."""
    pass


def _analyse(surface):
    src = np.asarray(surface.source_positions_mm, dtype=np.float64)
    pix = np.asarray(surface.pixel_centers_mm, dtype=np.float64)
    mask = np.asarray(surface.valid_mask)
    nv, nrows, ncols = mask.shape

    radius = np.hypot(src[:, 0], src[:, 1])
    if radius.min() <= 0 or np.ptp(radius) > 1e-4 * radius.mean() or np.ptp(src[:, 2]) > 1e-3:
        raise NotImplementedError(
            "fdk_surface needs a circular source orbit about z at constant height")

    col_ok = mask.all(axis=(0, 1))
    if not np.array_equal(mask.any(axis=(0, 1)), col_ok):
        raise NotImplementedError(
            "fdk_surface needs the same valid columns in every view and row (whole-column gaps only)")
    cols = np.flatnonzero(col_ok)
    if cols.size < 2:
        raise ValueError("fdk_surface needs at least two valid columns")

    p = pix[:, :, cols, :]
    h = p[..., :2] - src[:, None, None, :2]
    tol = 1e-3 * radius.mean() * 1e-2
    if np.ptp(h, axis=1).max() > tol or np.ptp(p[..., 2], axis=2).max() > tol \
            or np.ptp(p[..., 2], axis=0).max() > tol:
        raise NotImplementedError(
            "fdk_surface needs an extruded detector: pixel x/y independent of the row and "
            "z independent of the column and view")

    c = -src[:, :2] / radius[:, None]
    h0 = h[:, 0]
    gamma = np.arctan2(c[:, None, 0] * h0[..., 1] - c[:, None, 1] * h0[..., 0],
                       c[:, None, 0] * h0[..., 0] + c[:, None, 1] * h0[..., 1])
    rho = np.hypot(h0[..., 0], h0[..., 1])
    if np.ptp(gamma, axis=0).max() > 2e-5 or np.ptp(rho, axis=0).max() > 1e-4 * rho.mean():
        raise NotImplementedError("fdk_surface needs identical fan-angle tables in every view")
    gamma, rho = gamma[0], rho[0]
    zrow = p[0, :, 0, 2]

    d = np.diff(gamma)
    if np.all(d < 0):
        cols = cols[::-1]
        gamma, rho = gamma[::-1], rho[::-1]
    elif not np.all(d > 0):
        raise ValueError("fan angles must be strictly monotone along the columns")
    flip_rows = nrows > 1 and zrow[-1] < zrow[0]
    if flip_rows:
        zrow = zrow[::-1]
    if nrows > 1 and not np.all(np.diff(zrow) > 0):
        raise ValueError("detector rows must have strictly monotone z")

    fan = _Fan()
    fan.cols, fan.flip_rows = cols, flip_rows
    fan.gamma, fan.rho, fan.zrow = gamma, rho, zrow
    fan.radius, fan.sz, fan.src = radius.mean(), src[0, 2], src
    return fan


def _view_weights(beta):
    """Return (beta_relative >= 0 along the sweep, d_beta per view, sweep direction, span, is_full)."""
    beta = np.unwrap(beta)
    nv = beta.size
    if nv < 2:
        raise ValueError("fdk_surface needs at least two views")
    sweep = beta[-1] - beta[0]
    direction = 1.0 if sweep >= 0 else -1.0
    rel = direction * (beta - beta[0])
    steps = np.diff(rel)
    if np.any(steps <= 0):
        raise ValueError("the source must move monotonically along the orbit")
    span = rel[-1]
    tol = 1e-3
    if span > 2 * np.pi + tol:
        raise ValueError("fdk_surface supports at most one full orbit")
    dbeta = np.zeros(nv)
    if span >= 2 * np.pi - tol:
        dbeta[:-1] += steps / 2
        dbeta[1:] += steps / 2
        return rel, dbeta, direction, span, True
    if span + np.median(steps) >= 2 * np.pi - tol:
        # Periodic trapezoid rule, including the closing interval. A median
        # step loses the actual angular measure for nonuniform views.
        closing = 2 * np.pi - span
        dbeta[:-1] += steps / 2
        dbeta[1:] += steps / 2
        dbeta[[0, -1]] += closing / 2
        return rel, dbeta, direction, span, True
    dbeta[:-1] += steps / 2
    dbeta[1:] += steps / 2
    return rel, dbeta, direction, span, False
