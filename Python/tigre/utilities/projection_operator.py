"""Projection adapters for conventional TIGRE geometry and explicit surfaces.

Arrays use float32, volume order (z, y, x), and projection order
(views, rows, columns). ``subset`` indices refer to the current operator.
The conventional adapter uses TIGRE's legacy ``matched`` backprojector;
unlike the surface Siddon pair, it is not an exact discrete transpose.
"""

import copy
from typing import Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class ProjectionOperator(Protocol):
    """Solver contract; ``valid_mask=None`` means every sample is present.

    Implementations must not modify inputs and must return finite float32
    arrays of the declared shapes. Missing output samples must be zero;
    missing input samples must be ignored, including NaN/Inf there.
    Geometry, GPU selection and other backend state belong to the operator.
    Operators must be linear; row/column-normalized solvers additionally
    assume nonnegative projection weights, as in TIGRE's ray integrals.
    """

    volume_shape: tuple
    projection_shape: tuple
    valid_mask: object

    def forward(self, volume): ...
    def adjoint(self, projections): ...
    def subset(self, view_indices): ...


def _indices(indices, count):
    indices = np.asarray(indices)
    if indices.ndim != 1 or indices.size == 0 or indices.dtype.kind not in "iu":
        raise ValueError("view_indices must be a nonempty one-dimensional integer array")
    if np.any(indices < 0) or np.any(indices >= count):
        raise ValueError("view_indices are out of range")
    return indices.astype(np.intp)


def _array(value, shape, name, mask=None, *, copy_data=True):
    value = np.asarray(value)
    if np.iscomplexobj(value):
        raise TypeError(f"{name} must be real")
    if value.shape != tuple(shape):
        raise ValueError(f"{name} must have shape {tuple(shape)}, got {value.shape}")
    if copy_data or mask is not None:
        value = np.array(value, dtype=np.float32, order="C", copy=True)
    else:
        value = np.ascontiguousarray(value, dtype=np.float32)
    if mask is not None:
        value[~mask] = 0
    if not np.isfinite(value).all():
        raise ValueError(f"{name} contains NaN or Inf at valid samples")
    return value


def _volume_shape(geo):
    shape = np.asarray(geo.nVoxel)
    if shape.shape != (3,) or shape.dtype.kind not in "iu" or np.any(shape <= 0):
        raise ValueError("geo.nVoxel must contain three positive integers")
    return tuple(int(v) for v in shape)


def _device_indices(gpuids):
    """Validate explicit device selection without querying CUDA at binding time."""
    if gpuids is None:
        return None
    values = getattr(gpuids, "devices", gpuids)
    try:
        values = [values] if np.isscalar(values) else list(values)
        indices = np.asarray(values)
    except (TypeError, ValueError) as error:
        raise ValueError("gpuids must contain nonnegative integer device indices") from error
    if (
        indices.ndim != 1
        or not indices.size
        or indices.dtype.kind not in "iu"
        or np.any(indices < 0)
    ):
        raise ValueError("gpuids must contain nonnegative integer device indices")
    return [int(value) for value in indices]


def _conventional_gpuids(gpuids):
    if gpuids is None or hasattr(gpuids, "devices"):
        return copy.deepcopy(gpuids)
    from tigre.utilities.gpu import GpuIds

    devices = GpuIds()
    devices.devices = _device_indices(gpuids)
    return devices


def prepare_projections(projections, operator):
    """Own the data buffer and remove missing samples before arithmetic."""
    for name in ("volume_shape", "projection_shape"):
        try:
            raw = tuple(getattr(operator, name))
        except (AttributeError, TypeError) as error:
            raise TypeError(f"operator.{name} must contain three positive integers") from error
        if len(raw) != 3 or any(
            isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)) or v <= 0
            for v in raw
        ):
            raise ValueError(f"operator.{name} must contain three positive integers")
    if not all(
        callable(getattr(operator, name, None)) for name in ("forward", "adjoint", "subset")
    ):
        raise TypeError("operator must implement forward, adjoint and subset")
    mask = operator.valid_mask
    if mask is not None:
        mask = np.asarray(mask)
        if mask.dtype.kind != "b" or mask.shape != tuple(operator.projection_shape):
            raise ValueError("operator.valid_mask must be boolean with projection_shape")
    return _array(projections, operator.projection_shape, "projections", mask)


class TigreOperator:
    """Conventional Siddon / matched TIGRE pair with owned geometry and GPUs."""

    valid_mask = None

    def __init__(self, geo, angles, gpuids=None):
        _device_indices(gpuids)
        self.geo = copy.deepcopy(geo)
        angles = np.asarray(angles)
        if angles.size == 0 or not np.isfinite(angles).all():
            raise ValueError("angles must be nonempty and finite")
        self.geo.check_geo(angles)
        self.angles = self.geo.angles.copy()
        self.gpuids = copy.deepcopy(gpuids)
        self.volume_shape = _volume_shape(self.geo)
        self.projection_shape = (len(self.angles),) + tuple(int(v) for v in self.geo.nDetector)

    def forward(self, volume):
        from tigre.utilities.Ax import Ax

        return Ax(
            _array(volume, self.volume_shape, "volume", copy_data=False),
            self.geo,
            self.angles,
            "Siddon",
            gpuids=_conventional_gpuids(self.gpuids),
        )

    def adjoint(self, projections):
        from tigre.utilities.Atb import Atb

        return Atb(
            prepare_projections(projections, self),
            self.geo,
            self.angles,
            "matched",
            gpuids=_conventional_gpuids(self.gpuids),
        )

    def subset(self, view_indices):
        indices = _indices(view_indices, self.projection_shape[0])
        geo = copy.deepcopy(self.geo)
        # check_geo normalizes each field to a per-view array. Do not infer
        # fields from their length: nVoxel also has length three.
        for name in ("DSO", "DSD", "COR", "offOrigin", "offDetector", "rotDetector"):
            setattr(geo, name, getattr(geo, name)[indices].copy())
        return TigreOperator(geo, self.angles[indices], self.gpuids)


class SurfaceOperator:
    """Exact surface Siddon pair, including masks and optional GPU slabs."""

    def __init__(self, geo, surface, gpuids=None, n_slabs=None):
        _device_indices(gpuids)
        if n_slabs is not None and (
            isinstance(n_slabs, (bool, np.bool_))
            or not isinstance(n_slabs, (int, np.integer))
            or n_slabs < 1
        ):
            raise ValueError("n_slabs must be a positive integer or None")
        self.geo = copy.deepcopy(geo)
        self.surface = copy.deepcopy(surface)
        self.gpuids = copy.deepcopy(gpuids)
        self.n_slabs = n_slabs
        self.volume_shape = _volume_shape(self.geo)
        # Reject view-dependent origins before taking a subset could hide
        # their variation. Surface kernels require a single volume frame.
        offset = np.asarray(self.geo.offOrigin)
        if offset.ndim == 2 and not np.all(offset == offset[0]):
            raise ValueError("offOrigin must be identical for all views")
        voxel = np.asarray(self.geo.dVoxel)
        if (
            voxel.shape != (3,)
            or not np.isfinite(voxel).all()
            or np.any(voxel <= 0)
            or offset.shape not in ((3,), (1, 3), (surface.valid_mask.shape[0], 3))
            or not np.isfinite(offset).all()
        ):
            raise ValueError("surface volume spacing/origin must be finite with positive dVoxel")
        self.projection_shape = tuple(self.surface.valid_mask.shape)
        self.valid_mask = self.surface.valid_mask

    def forward(self, volume):
        from tigre.utilities.Ax_surface import Ax_surface

        return Ax_surface(
            _array(volume, self.volume_shape, "volume", copy_data=False),
            self.geo,
            self.surface,
            gpuids=copy.deepcopy(self.gpuids),
            n_slabs=self.n_slabs,
        )

    def adjoint(self, projections):
        from tigre.utilities.Ax_surface import Atb_surface

        return Atb_surface(
            prepare_projections(projections, self),
            self.geo,
            self.surface,
            gpuids=copy.deepcopy(self.gpuids),
            n_slabs=self.n_slabs,
        )

    def subset(self, view_indices):
        indices = _indices(view_indices, self.projection_shape[0])
        geo = copy.deepcopy(self.geo)
        if np.asarray(geo.offOrigin).ndim == 2:
            geo.offOrigin = np.asarray(geo.offOrigin)[0].copy()
        # DetectorSurfaceGeometry.subset already owns its arrays. Avoid a
        # second acquisition-sized deep copy during every subset construction.
        result = object.__new__(type(self))
        result.geo = geo
        result.surface = self.surface.subset(indices)
        result.gpuids = copy.deepcopy(self.gpuids)
        result.n_slabs = self.n_slabs
        result.volume_shape = self.volume_shape
        result.projection_shape = tuple(result.surface.valid_mask.shape)
        result.valid_mask = result.surface.valid_mask
        return result
