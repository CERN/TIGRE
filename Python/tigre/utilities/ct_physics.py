"""Physical transmission-CT checks, separate from projector/FDK restrictions."""
import numpy as np


def _segment_hits_box(source, pixels, lower, upper):
    """Positive-length intersections of finite source-to-pixel segments."""
    enter = np.zeros(len(pixels))
    leave = np.ones(len(pixels))
    for axis in range(3):
        direction = pixels[:, axis] - source[axis]
        parallel = direction == 0
        a = np.divide(lower[axis] - source[axis], direction,
                      out=np.full(len(pixels), -np.inf), where=~parallel)
        b = np.divide(upper[axis] - source[axis], direction,
                      out=np.full(len(pixels), np.inf), where=~parallel)
        near, far = np.minimum(a, b), np.maximum(a, b)
        if source[axis] < lower[axis] or source[axis] > upper[axis]:
            near[parallel], far[parallel] = np.inf, -np.inf
        enter = np.maximum(enter, near)
        leave = np.minimum(leave, far)
    return leave > enter + 1e-12


def validate_transmission_geometry(surface, geometry):
    """Check hardware positions and require some transmission over the scan.

    Source and valid detector pixels must lie outside the volume. Empty views
    and panels are allowed; only an entirely nontransmitting scan is rejected.
    This accommodates helical entry/exit and overscan. This does not
    certify detector collision clearance, sampling completeness or FDK validity.
    """
    if geometry is None:
        raise ValueError("Volume size and center are required to check physical CT geometry.")
    if hasattr(geometry, "object_angles_deg"):
        from .object_motion import fixed_object_geometry
        surface, geometry = fixed_object_geometry(surface, geometry)
    from .visualization.surface_renderer import _volume_box
    _volume_box(geometry)
    source = np.asarray(surface.source_positions_mm)
    centers = np.atleast_2d(np.asarray(geometry.offOrigin, dtype=float))[:, ::-1]
    if len(centers) not in (1, len(source)):
        raise ValueError("offOrigin must contain one offset or one per view.")
    half = np.asarray(geometry.sVoxel, dtype=float)[::-1] / 2
    tolerance = 1e-9 * max(float(half.max()), 1)
    modules = getattr(surface, "_module_slices", None) or [
        (slice(None), slice(None))]
    measured, empty_views = 0, 0
    for view, position in enumerate(source):
        center = centers[0 if len(centers) == 1 else view]
        lower, upper = center - half, center + half
        relative_source = position - center
        if np.all(np.abs(relative_source) <= half + tolerance):
            raise ValueError(f"Physical CT geometry: source is inside or on the volume at view {view + 1}.")
        view_hits = 0
        for panel, (rows, columns) in enumerate(modules, 1):
            grid = surface.pixel_centers_mm[view, rows, columns]
            valid = surface.valid_mask[view, rows, columns].ravel()
            pixels = grid.reshape(-1, 3)
            panel_hits = 0
            for start in range(0, len(pixels), 16384):
                values = pixels[start:start + 16384][valid[start:start + 16384]]
                if not len(values):
                    continue
                relative = values - center
                if np.any(np.all(np.abs(relative) <= half + tolerance, axis=1)):
                    raise ValueError(f"Physical CT geometry: detector panel {panel} is inside or on the volume at view {view + 1}.")
                panel_hits += int(np.count_nonzero(_segment_hits_box(position, values, lower, upper)))
            view_hits += panel_hits
        if not view_hits:
            empty_views += 1
        measured += view_hits
    if not measured:
        raise ValueError(
            "Physical CT geometry: no measured source-to-detector ray crosses the volume "
            "in the entire scan. Check source/detector sides, volume position, scan span and pixel sampling.")
    return {"views": len(source), "crossing_samples": measured, "empty_views": empty_views}
