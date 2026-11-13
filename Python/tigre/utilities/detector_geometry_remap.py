"""Projection remapping for approximate surface FDK reconstruction."""

import numpy as np

from .detector_geometry import DetectorSurfaceGeometry


def remap_projections_for_fdk(projections, surface, target_columns=None,
                              target_row_positions_mm=None):
    """Resample projections onto a circular, equiangular extruded detector.

    Input projections must have shape ``(views, rows, columns)`` and match the
    supplied surface. Interpolation is performed in fan-angle and vertical
    ray-angle coordinates, separately within contiguous valid column segments.
    Samples outside measured coverage remain masked; no extrapolation is done.
    target_row_positions_mm optionally specifies virtual detector z coordinates
    in world mm. Choose an interior range for tilted panels to avoid requiring
    unsupported edge rays across every row. This does not fill missing coverage.

    Returns ``(remapped_projections, target_surface, report)``. This is an
    approximation and cannot make a non-circular source trajectory suitable
    for FDK.
    """
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import QhullError

    source = np.asarray(surface.source_positions_mm, dtype=np.float64)
    pixels = np.asarray(surface.pixel_centers_mm, dtype=np.float64)
    valid = np.asarray(surface.valid_mask, dtype=bool)
    data = np.asarray(projections, dtype=np.float32)
    if data.shape != valid.shape:
        raise ValueError(
            f"projections must have shape {valid.shape}; got {data.shape}")
    if not np.isfinite(data[valid]).all():
        raise ValueError("projections contain NaN or Inf at valid detector samples")

    radius = np.hypot(source[:, 0], source[:, 1])
    if (radius.min() <= 0 or np.ptp(radius) > 1e-4 * radius.mean()
            or np.ptp(source[:, 2]) > 1e-3):
        raise ValueError(
            "Remapping for fdk_surface requires a circular source orbit "
            "about z at constant height; detector interpolation cannot repair "
            "a non-circular source path.")
    beta = np.unwrap(np.arctan2(source[:, 1], source[:, 0]))
    if np.any(np.diff(beta) == 0) or not (
            np.all(np.diff(beta) > 0) or np.all(np.diff(beta) < 0)):
        raise ValueError("Source positions must move monotonically around the orbit.")

    nv, nrows, ncols = valid.shape
    if not np.array_equal(valid, np.broadcast_to(valid[0], valid.shape)):
        raise ValueError(
            "Remapping requires the same detector mask in every view.")
    whole_columns = valid.all(axis=1)
    if not np.array_equal(whole_columns, valid.any(axis=1)):
        raise ValueError(
            "Remapping requires whole-column detector gaps; partial-column "
            "pixel masks are not supported.")
    valid_columns = whole_columns[0]
    segments = []
    column_indices = np.flatnonzero(valid_columns)
    if column_indices.size < 2:
        raise ValueError("At least two valid detector columns are required.")
    module_slices = getattr(surface, "_module_slices", None)
    if module_slices:
        segments = []
        for rows, columns in module_slices:
            if rows.start != 0 or rows.stop != nrows:
                raise ValueError(
                    "Remapping requires detector modules to span the same rows.")
            segments.append(np.arange(columns.start, columns.stop))
        covered = np.zeros(ncols, dtype=bool)
        for part in segments:
            if covered[part].any() or not valid_columns[part].all():
                raise ValueError(
                    "Module slices overlap or include invalid detector columns.")
            covered[part] = True
        if not np.array_equal(covered, valid_columns):
            raise ValueError("Module slices do not match the valid detector columns.")
    else:
        breaks = np.flatnonzero(np.diff(column_indices) > 1) + 1
        segments = [part for part in np.split(column_indices, breaks) if part.size]
    if any(part.size < 2 for part in segments):
        raise ValueError(
            "Each contiguous detector segment needs at least two columns.")

    def ray_angles(view):
        inward = -source[view, :2] / radius[view]
        delta = pixels[view, ..., :2] - source[view, None, None, :2]
        fan = np.arctan2(
            inward[0] * delta[..., 1] - inward[1] * delta[..., 0],
            inward[0] * delta[..., 0] + inward[1] * delta[..., 1],
        )
        radial = np.hypot(delta[..., 0], delta[..., 1])
        vertical = np.arctan2(
            pixels[view, ..., 2] - source[view, None, None, 2], radial)
        if not np.isfinite(fan).all() or not np.isfinite(vertical).all():
            raise ValueError("Detector rays must have finite fan and vertical angles.")
        return fan, radial, vertical

    # Keep only one view's ray-coordinate tables in memory at a time. This is
    # important for large acquisitions such as the 5,701-view GCT scan.
    lower_by_view = np.empty(nv)
    upper_by_view = np.empty(nv)
    pitch_by_view = []
    radius_by_view = np.empty(nv)
    row_z_by_view = np.empty((nv, nrows))
    for view in range(nv):
        fan, radial, _ = ray_angles(view)
        lower_by_view[view] = fan[:, valid_columns].min()
        upper_by_view[view] = fan[:, valid_columns].max()
        radius_by_view[view] = np.median(radial[:, valid_columns])
        row_z_by_view[view] = np.nanmedian(
            np.where(valid_columns[None, :], pixels[view, ..., 2], np.nan),
            axis=1)
        for indices in segments:
            for row in range(nrows):
                differences = np.diff(fan[row, indices])
                if not (np.all(differences > 0) or np.all(differences < 0)):
                    raise ValueError(
                        "Each detector segment must have strictly monotone "
                        "fan angles along its columns.")
                pitch_by_view.append(float(np.median(np.abs(differences))))
    angular_pitch = float(np.median(pitch_by_view))
    if not np.isfinite(angular_pitch) or angular_pitch <= 0:
        raise ValueError("Could not determine a positive detector fan-angle pitch.")

    lower = float(lower_by_view.max())
    upper = float(upper_by_view.min())
    if not upper > lower:
        raise ValueError("Detector views have no common fan-angle coverage.")
    if target_columns is None:
        target_columns = int(np.floor((upper - lower) / angular_pitch)) + 1
    if (not isinstance(target_columns, (int, np.integer))
            or isinstance(target_columns, (bool, np.bool_))
            or target_columns < 2):
        raise ValueError("target_columns must be an integer of at least two.")
    if target_row_positions_mm is None:
        target_z = np.median(row_z_by_view, axis=0)
    else:
        target_z = np.asarray(target_row_positions_mm, dtype=np.float64)
        if target_z.ndim != 1 or target_z.size < 1:
            raise ValueError("target_row_positions_mm must be a nonempty 1D array.")
    target_nrows = target_z.size
    if nrows == 1 and target_nrows != 1:
        raise ValueError("A one-row input requires a one-row target.")
    if nv * target_nrows * int(target_columns) > 20_000_000:
        raise ValueError(
            "Remapped data would exceed the 20,000,000-sample dense surface limit; "
            "reduce the input views or detector size.")
    target_gamma = np.linspace(lower, upper, int(target_columns))

    target_radius = float(np.median(radius_by_view))
    if not np.isfinite(target_radius) or target_radius <= 0:
        raise ValueError("Could not determine a positive target detector radius.")
    if not np.isfinite(target_z).all():
        raise ValueError("Could not determine target detector row positions.")
    if target_nrows > 1 and not (
            np.all(np.diff(target_z) > 0) or np.all(np.diff(target_z) < 0)):
        raise ValueError("Target detector rows must have strictly monotone z positions.")

    # Interpolate only inside each contiguous measured segment. This prevents
    # a triangulation from silently filling physical gaps between panels.
    remapped = np.zeros((nv, target_nrows, int(target_columns)), dtype=np.float32)
    supported = np.zeros(remapped.shape, dtype=bool)
    for view in range(nv):
        fan, _, vertical = ray_angles(view)
        query_gamma, query_eta = np.meshgrid(
            target_gamma,
            np.arctan2(target_z - source[view, 2], target_radius),
        )
        overlap_count = np.zeros(remapped.shape[1:], dtype=np.uint8)
        for indices in segments:
            points_fan = fan[:, indices]
            points_vertical = vertical[:, indices]
            samples = data[view][:, indices]
            query_min = float(points_fan.min())
            query_max = float(points_fan.max())
            first = int(np.searchsorted(target_gamma, query_min, side="left"))
            last = int(np.searchsorted(target_gamma, query_max, side="right"))
            columns = np.arange(first, last)
            if columns.size < 1:
                continue

            if nrows == 1:
                eta = points_vertical[0]
                if np.ptp(eta) > 1e-6:
                    raise ValueError(
                        "A one-row detector segment has varying vertical ray "
                        "angles; safe 1D remapping is not possible.")
                if abs(float(np.median(query_eta[0, columns])) - float(np.median(eta))) > 1e-5:
                    raise ValueError(
                        "Target fan rays do not match the one-row detector's "
                        "vertical ray angle.")
                order = np.argsort(points_fan[0])
                values = np.interp(
                    target_gamma[columns], points_fan[0, order], samples[0, order],
                    left=np.nan, right=np.nan)
                remapped[view, 0, columns] = np.nan_to_num(values)
                supported[view, 0, columns] |= np.isfinite(values)
                overlap_count[0, columns] += np.isfinite(values)
                continue

            points = np.column_stack((
                points_fan.ravel(), points_vertical.ravel()))
            queries = np.column_stack((
                query_gamma[:, columns].ravel(),
                query_eta[:, columns].ravel()))
            try:
                interpolator = LinearNDInterpolator(points, samples.ravel(), fill_value=np.nan)
            except (QhullError, ValueError) as exc:
                raise ValueError(
                    "Could not triangulate a detector segment for interpolation: "
                    f"{exc}") from exc
            values = np.asarray(interpolator(queries)).reshape(target_nrows, -1)
            finite = np.isfinite(values)
            remapped[view][:, columns] = np.where(
                finite, values, remapped[view][:, columns])
            supported[view][:, columns] |= finite
            overlap_count[:, columns] += finite
        if np.any(overlap_count > 1):
            raise ValueError(
                "Multiple detector segments cover the same target rays; "
                "automatic remapping is ambiguous for this layout.")

    common_columns = supported.all(axis=(0, 1))
    if common_columns.sum() < 2:
        raise ValueError(
            "Fewer than two target columns are measured in every view and row; "
            "there is insufficient common detector coverage.")
    target_mask = np.broadcast_to(
        common_columns[None, None, :], remapped.shape).copy()
    target_source = source.astype(np.float32)
    target_pixels = np.empty((nv, target_nrows, int(target_columns), 3), dtype=np.float32)
    for view in range(nv):
        inward_xy = -source[view, :2] / radius[view]
        tangent_xy = np.array([-inward_xy[1], inward_xy[0]])
        direction_xy = (
            np.cos(target_gamma)[:, None] * inward_xy
            + np.sin(target_gamma)[:, None] * tangent_xy)
        target_pixels[view, :, :, :2] = (
            source[view, None, None, :2]
            + target_radius * direction_xy[None, :, :])
        target_pixels[view, :, :, 2] = target_z[:, None]

    target_surface = DetectorSurfaceGeometry(
        target_source, target_pixels, target_mask)
    from tigre.utilities._surface_geometry_checks import _analyse, _view_weights

    target_fan = _analyse(target_surface)
    rel, _, _, span, full_scan = _view_weights(np.arctan2(
        target_source[:, 1], target_source[:, 0]))
    del rel
    if full_scan:
        short_scan_deficiency = 0.0
    else:
        required_gamma = 0.5 * (span - np.pi)
        short_scan_deficiency = max(
            0.0, float(np.abs(target_fan.gamma).max() - required_gamma))
    report = {
        "input_shape": tuple(data.shape),
        "output_shape": tuple(remapped.shape),
        "target_radius_mm": target_radius,
        "target_row_positions_mm": target_z.tolist(),
        "panel_count": len(segments),
        "supported_column_count": int(common_columns.sum()),
        "missing_column_count": int((~common_columns).sum()),
        "fan_angle_range_rad": (float(target_gamma[0]), float(target_gamma[-1])),
        "short_scan_supported": bool(
            full_scan or short_scan_deficiency <= np.deg2rad(2.0)),
        "short_scan_deficiency_deg": float(np.degrees(short_scan_deficiency)),
        "interpolation": "linear in fan and vertical ray angles; no extrapolation",
        "warning": (
            "Projection values were interpolated onto an idealized circular "
            "equiangular detector. This changes measured samples and may blur "
            "or bias data; validate reconstruction quality against the original "
            "measurements."),
    }
    remapped[~target_mask] = 0.0
    return remapped, target_surface, report
