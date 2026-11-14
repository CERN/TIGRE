"""Explicit source/pixel geometry and rigid detector modules. Coordinates are world (x, y, z) in mm."""

import numpy as np

__all__ = ["DetectorModule", "DetectorSurfaceGeometry"]


def _source_positions(source_positions_mm):
    source = np.asarray(source_positions_mm, dtype=np.float64)
    if source.ndim != 2 or source.shape[1] != 3 or source.shape[0] == 0:
        raise ValueError("source_positions_mm must have shape (n_views, 3), n_views > 0")
    if not np.isfinite(source).all():
        raise ValueError("source_positions_mm must contain only finite values")
    source = np.ascontiguousarray(source, dtype=np.float32)
    if not np.isfinite(source).all():
        raise ValueError("source_positions_mm is out of float32 range")
    return source


def _pixel_centers(pixel_centers_mm, n_views):
    centers = np.asarray(pixel_centers_mm)
    if centers.dtype != np.float32:
        centers = np.asarray(centers, dtype=np.float64)
    if centers.ndim == 3 and centers.shape[-1] == 3:
        centers = np.broadcast_to(centers, (n_views,) + centers.shape)
    if (
        centers.ndim != 4
        or centers.shape[0] != n_views
        or centers.shape[-1] != 3
        or (min(centers.shape[1:3]) <= 0)
    ):
        raise ValueError(
            "pixel_centers_mm must have shape (rows, columns, 3) or (n_views, rows, columns, 3)"
        )
    if not np.isfinite(centers).all():
        raise ValueError("pixel_centers_mm must contain only finite values")
    centers = np.ascontiguousarray(centers, dtype=np.float32)
    if not np.isfinite(centers).all():
        raise ValueError("pixel_centers_mm is out of float32 range")
    return centers


def _pixel_mask(mask, shape):
    n_views, n_rows, n_columns = shape
    if mask is None:
        return np.ones(shape, dtype=np.bool_)
    values = np.asarray(mask)
    if values.dtype.kind not in "biuf" or np.any((values != 0) & (values != 1)):
        raise ValueError("valid_mask must contain only boolean or 0/1 values")
    valid = values.astype(np.bool_, copy=False)
    if valid.shape == (n_rows, n_columns):
        valid = np.broadcast_to(valid, shape)
    if valid.shape != shape:
        raise ValueError(
            f"valid_mask must have shape ({n_rows}, {n_columns}) or ({n_views}, {n_rows}, {n_columns}); got {valid.shape}"
        )
    return np.ascontiguousarray(valid)


class DetectorModule:
    """A detector module with local pixel centers and per-view rigid poses.

    ``pixel_centers_mm`` is shaped ``(rows, columns, 3)`` and may describe
    either a flat or curved panel, including nonuniform pixel spacing.
    ``pose_matrices`` maps local homogeneous coordinates to the lab frame and
    may be one ``(4, 4)`` transform or one transform per view. Integer offsets
    place the module in the common dense detector-array layout; unoccupied
    cells become invalid gap samples.
    """

    def __init__(self, pixel_centers_mm, pose_matrices, row_offset=0, column_offset=0):
        centers = np.asarray(pixel_centers_mm, dtype=np.float64)
        if centers.ndim != 3 or centers.shape[-1] != 3 or min(centers.shape[:2]) <= 0:
            raise ValueError("module pixel_centers_mm must have shape (rows, columns, 3)")
        if not np.isfinite(centers).all():
            raise ValueError("module pixel_centers_mm must contain only finite values")
        if (
            not isinstance(row_offset, (int, np.integer))
            or isinstance(row_offset, (bool, np.bool_))
            or row_offset < 0
        ):
            raise ValueError("row_offset must be a non-negative integer")
        if (
            not isinstance(column_offset, (int, np.integer))
            or isinstance(column_offset, (bool, np.bool_))
            or column_offset < 0
        ):
            raise ValueError("column_offset must be a non-negative integer")
        transforms = np.asarray(pose_matrices, dtype=np.float64)
        if transforms.shape == (4, 4):
            transforms = transforms[None, ...]
        if transforms.ndim != 3 or transforms.shape[1:] != (4, 4):
            raise ValueError("pose_matrices must have shape (4, 4) or (n_views, 4, 4)")
        if not np.isfinite(transforms).all():
            raise ValueError("pose_matrices must contain only finite values")
        if not np.allclose(transforms[:, 3, :], [0.0, 0.0, 0.0, 1.0], atol=1e-07):
            raise ValueError("pose_matrices must be homogeneous rigid transforms")
        rotations = transforms[:, :3, :3]
        identity = np.eye(3)
        if not np.allclose(
            rotations @ np.swapaxes(rotations, 1, 2), identity, atol=1e-05
        ) or np.any(np.linalg.det(rotations) <= 0):
            raise ValueError("pose_matrices must contain right-handed rigid rotations")
        self.pixel_centers_mm = np.ascontiguousarray(centers, dtype=np.float32)
        self.pose_matrices = np.ascontiguousarray(transforms, dtype=np.float32)
        if (
            not np.isfinite(self.pixel_centers_mm).all()
            or not np.isfinite(self.pose_matrices).all()
        ):
            raise ValueError("module coordinates or poses exceed float32 range")
        self.row_offset = int(row_offset)
        self.column_offset = int(column_offset)

    def world_pixel_centers(self, n_views):
        if (
            not isinstance(n_views, (int, np.integer))
            or isinstance(n_views, (bool, np.bool_))
            or n_views < 1
        ):
            raise ValueError("n_views must be a positive integer")
        if self.pose_matrices.shape[0] not in (1, n_views):
            raise ValueError("module pose count must be one or match the source view count")
        poses = np.broadcast_to(self.pose_matrices, (n_views, 4, 4))
        local = np.concatenate(
            (
                self.pixel_centers_mm,
                np.ones(self.pixel_centers_mm.shape[:2] + (1,), dtype=np.float32),
            ),
            axis=2,
        )
        return np.einsum("vij,hwj->vhwi", poses, local)[..., :3]


class DetectorSurfaceGeometry:
    """Per-view source positions, detector pixel centres and a validity mask.

    Geometry coordinates are world (x, y, z) in millimetres. Pixel arrays use
    (views, rows, columns, 3); projection arrays use (views, rows, columns).
    False mask entries are absent measurements, not measured zero values.
    This model does not implement projection, reconstruction or gap filling.
    """

    def __init__(self, source_positions_mm, pixel_centers_mm, valid_mask=None):
        self.source_positions_mm = _source_positions(source_positions_mm)
        self.pixel_centers_mm = _pixel_centers(pixel_centers_mm, len(self.source_positions_mm))
        self.valid_mask = _pixel_mask(valid_mask, self.pixel_centers_mm.shape[:3])
        self._module_slices = None

    def subset(self, view_indices):
        """Return a surface restricted to the given views (data is copied)."""
        idx = np.asarray(view_indices)
        if idx.ndim != 1 or idx.size == 0 or idx.dtype.kind not in "iu":
            raise ValueError("view_indices must be a nonempty integer vector")
        if np.any(idx < 0) or np.any(idx >= len(self.source_positions_mm)):
            raise ValueError("view_indices are out of range")
        idx = idx.astype(np.intp)
        out = object.__new__(DetectorSurfaceGeometry)
        out.source_positions_mm = np.ascontiguousarray(self.source_positions_mm[idx])
        out.pixel_centers_mm = np.ascontiguousarray(self.pixel_centers_mm[idx])
        out.valid_mask = np.ascontiguousarray(self.valid_mask[idx])
        out._module_slices = self._module_slices
        return out

    @classmethod
    def from_modules(cls, source_positions_mm, modules):
        """Combine modules into one dense array, leaving gaps explicitly masked."""
        source = _source_positions(source_positions_mm)
        modules = tuple(modules)
        if not modules or not all((isinstance(m, DetectorModule) for m in modules)):
            raise ValueError("modules must be a non-empty sequence of DetectorModule objects")
        n_rows = 0
        n_columns = 0
        for module in modules:
            row_end = module.row_offset + module.pixel_centers_mm.shape[0]
            column_end = module.column_offset + module.pixel_centers_mm.shape[1]
            n_rows = max(n_rows, row_end)
            n_columns = max(n_columns, column_end)
        owner = np.zeros((n_rows, n_columns), dtype=np.bool_)
        module_slices = []
        for module in modules:
            row_end = module.row_offset + module.pixel_centers_mm.shape[0]
            column_end = module.column_offset + module.pixel_centers_mm.shape[1]
            rows = slice(module.row_offset, row_end)
            columns = slice(module.column_offset, column_end)
            if owner[rows, columns].any():
                raise ValueError("detector modules overlap in the dense data layout")
            owner[rows, columns] = True
            module_slices.append((rows, columns))
        dense_centers = np.zeros((len(source), n_rows, n_columns, 3), dtype=np.float32)
        dense_mask = np.zeros((n_rows, n_columns), dtype=np.bool_)
        for module, (rows, columns) in zip(modules, module_slices):
            dense_centers[:, rows, columns, :] = module.world_pixel_centers(len(source))
            dense_mask[rows, columns] = True
        result = cls(source, dense_centers, dense_mask)
        result._module_slices = tuple(module_slices)
        return result

    def pack_module_data(self, module_data, fill_value=0.0):
        """Pack per-module ``(views, rows, columns)`` values into the dense layout."""
        if self._module_slices is None:
            raise ValueError("pack_module_data is available only for module geometry")
        if len(module_data) != len(self._module_slices):
            raise ValueError("module_data count must match the module geometry")
        try:
            fill = float(fill_value)
        except (TypeError, ValueError) as exc:
            raise ValueError("fill_value must be a finite scalar") from exc
        if not np.isfinite(fill):
            raise ValueError("fill_value must be a finite scalar")
        packed = np.full(self.valid_mask.shape, fill, dtype=np.float32)
        for data, (rows, columns) in zip(module_data, self._module_slices):
            values = np.asarray(data)
            expected = (
                len(self.source_positions_mm),
                rows.stop - rows.start,
                columns.stop - columns.start,
            )
            if values.shape != expected:
                raise ValueError(f"module data must have shape {expected}; got {values.shape}")
            if not np.isfinite(values).all():
                raise ValueError("module data must contain only finite values")
            packed[:, rows, columns] = values
        return packed

    def unpack_module_data(self, dense_data):
        """Extract a dense ``(views, rows, columns)`` array into module arrays."""
        if self._module_slices is None:
            raise ValueError("unpack_module_data is available only for module geometry")
        values = np.asarray(dense_data)
        if values.shape != self.valid_mask.shape:
            raise ValueError(
                f"dense_data must have shape {self.valid_mask.shape}; got {values.shape}"
            )
        if not np.isfinite(values).all():
            raise ValueError("dense_data must contain only finite values")
        return [
            np.ascontiguousarray(values[:, rows, columns]) for rows, columns in self._module_slices
        ]
