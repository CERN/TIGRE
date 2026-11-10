import numpy as np
from _AxSurface import _Ax_surface_ext, _Atb_surface_ext
from .gpu import GpuIds


def _volume_params(geo):
    counts = np.asarray(geo.nVoxel)
    if (
        counts.shape != (3,)
        or counts.dtype.kind not in "iu"
        or np.any(counts < 1)
        or np.any(counts > np.iinfo(np.int32).max)
    ):
        raise ValueError("geo.nVoxel must contain three positive int32 dimensions")
    offset = np.asarray(geo.offOrigin, dtype=np.float64)
    if offset.ndim == 2:
        if offset.shape[0] == 0 or not np.all(offset == offset[0]):
            raise ValueError("offOrigin must be identical for all views")
        offset = offset[0]
    if offset.shape != (3,) or not np.isfinite(offset).all():
        raise ValueError("offOrigin must contain three finite values")
    voxel = np.asarray(geo.dVoxel, dtype=np.float64)
    if voxel.shape != (3,) or not np.isfinite(voxel).all() or np.any(voxel <= 0):
        raise ValueError("dVoxel must contain three positive finite values")
    return (voxel[2], voxel[1], voxel[0]), (offset[2], offset[1], offset[0])


def _surface_args(surface):
    mask = surface.valid_mask
    if mask.all():
        mask = None
    return surface.source_positions_mm, surface.pixel_centers_mm, mask


def _devices(gpuids):
    """Device index list from a GpuIds, an int, an iterable of ints, or None (all GPUs)."""
    if gpuids is None:
        gpuids = GpuIds()
    if isinstance(gpuids, GpuIds):
        values = gpuids.devices
    elif np.ndim(gpuids) == 0:
        values = [gpuids]
    else:
        values = gpuids
    values = np.asarray(values)
    if (
        values.ndim != 1
        or values.size == 0
        or values.dtype.kind not in "iu"
        or np.any(values < 0)
        or np.any(values > np.iinfo(np.int32).max)
    ):
        raise ValueError("gpuids must contain nonnegative int32 device indices")
    return [int(d) for d in values]


def _slabs(nz, devices, n_slabs):
    if n_slabs is not None and (
        not isinstance(n_slabs, (int, np.integer))
        or isinstance(n_slabs, (bool, np.bool_))
        or n_slabs < 1
    ):
        raise ValueError("n_slabs must be a positive integer")
    n = len(devices) if n_slabs is None else int(n_slabs)
    n = min(n, nz)
    edges = np.linspace(0, nz, n + 1).astype(int)
    return [(int(edges[i]), int(edges[i + 1]), devices[i % len(devices)]) for i in range(n)]


def _slab_offset(offset, voxel, nz, lo, hi):
    # z is the last entry of the (x, y, z) tuples; slab centre shift in mm
    shift = (lo + hi - nz) * 0.5 * voxel[2]
    return (offset[0], offset[1], offset[2] + shift)


def _run_slabs(slabs, func):
    if len(slabs) == 1:
        return [func(slabs[0])]
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=len(slabs)) as pool:
        return list(pool.map(func, slabs))


def Atb_surface(proj, geo, surface, gpuids=None, n_slabs=None):
    """Exact adjoint of :func:`Ax_surface` (same Siddon ray weights).

    ``proj`` is float32 ``(views, rows, columns)``; the result has shape
    ``geo.nVoxel``. Masked samples contribute nothing. Float atomics make the
    result reproducible only up to summation-order rounding.

    ``gpuids`` is a :class:`tigre.utilities.gpu.GpuIds` (default: all GPUs, as in
    ``tigre.Atb``). The volume is split into z-slabs (``n_slabs``, default one per
    device) processed concurrently, so it need not fit on one GPU.
    """
    if proj.dtype != np.float32:
        raise TypeError("Input data should be float32, not " + str(proj.dtype))
    source, pixels, mask = _surface_args(surface)
    expected = (len(source),) + tuple(pixels.shape[1:3])
    if tuple(proj.shape) != expected:
        raise ValueError(
            "proj must have shape (views, rows, columns) = "
            + str(expected)
            + ", not "
            + str(proj.shape)
        )
    voxel, offset = _volume_params(geo)
    nvox = tuple(int(v) for v in geo.nVoxel)
    proj = np.array(proj, dtype=np.float32, order="C", copy=True)
    proj[~surface.valid_mask] = 0
    if not np.isfinite(proj).all():
        raise ValueError("valid projection samples must be finite")
    slabs = _slabs(nvox[0], _devices(gpuids), n_slabs)

    def run(slab):
        lo, hi, dev = slab
        return _Atb_surface_ext(
            proj,
            (hi - lo,) + nvox[1:],
            voxel,
            _slab_offset(offset, voxel, nvox[0], lo, hi),
            source,
            pixels,
            mask,
            dev,
        )

    parts = _run_slabs(slabs, run)
    return parts[0] if len(parts) == 1 else np.concatenate(parts, axis=0)


def Ax_surface(img, geo, surface, gpuids=None, n_slabs=None):
    """Forward projection onto an arbitrary detector surface (Siddon, cone beam).

    Parameters
    ----------
    img : float32 array of shape ``geo.nVoxel`` (z, y, x).
    geo : TIGRE geometry; only the volume fields (``nVoxel``, ``dVoxel``,
        ``offOrigin``) are used. ``offOrigin`` must be the same for all views.
    surface : DetectorSurfaceGeometry with source positions and pixel centers
        (mm) expressed in the volume-fixed lab frame.
    gpuids : :class:`tigre.utilities.gpu.GpuIds`, default all GPUs (as in
        ``tigre.Ax``). The volume is split into z-slabs run concurrently and
        the line integrals are summed, so it need not fit on one GPU.

    Returns
    -------
    float32 array ``(views, rows, columns)``; masked samples are 0.
    """
    if img.dtype != np.float32:
        raise TypeError("Input data should be float32, not " + str(img.dtype))
    if tuple(img.shape) != tuple(np.asarray(geo.nVoxel, dtype=int)):
        raise ValueError(
            "Input data should be of shape geo.nVoxel: "
            + str(geo.nVoxel)
            + " not:"
            + str(img.shape)
        )
    voxel, offset = _volume_params(geo)
    source, pixels, mask = _surface_args(surface)
    nz = img.shape[0]
    slabs = _slabs(nz, _devices(gpuids), n_slabs)

    def run(slab):
        lo, hi, dev = slab
        return _Ax_surface_ext(
            np.ascontiguousarray(img[lo:hi]),
            voxel,
            _slab_offset(offset, voxel, nz, lo, hi),
            source,
            pixels,
            mask,
            dev,
        )

    parts = _run_slabs(slabs, run)
    out = parts[0]
    for p in parts[1:]:
        out += p
    return out
