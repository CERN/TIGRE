"""Headless modular detector recipes, timing and compatibility checks.

Kept in TIGRE so geometry archives and numerical workflows do not depend on
the separately installed tigre-geometry-design-tool desktop application.
"""

import numpy as np
from tigre.utilities.detector_geometry import (
    cylindrical_modular_detector_surface,
    linear_source_surface,
    piecewise_linear_modular_detector_surface,
)

MAX_RECIPE_VIEWS = 6000

MAX_MODULE_DIMENSION = 4096

MAX_MODULES = 256

MAX_TOTAL_PIXELS = 20_000_000


def acquisition_timing(n_views, acquisition):
    """One synchronized detector frame per source pulse, starting at t=0."""
    names = (
        "readout_period_ms",
        "firing_rate_hz",
        "rotation_period_s",
        "table_speed_mm_s",
        "playback_rate",
    )
    if not isinstance(acquisition, dict) or set(acquisition) != set(names):
        raise ValueError(
            "Acquisition requires readout period, firing rate, rotation period, table speed and playback rate."
        )
    values = {name: float(acquisition[name]) for name in names}
    if not np.isfinite(list(values.values())).all():
        raise ValueError("Acquisition parameters must be finite.")
    if any(values[name] <= 0 for name in names if name != "table_speed_mm_s"):
        raise ValueError("Timing periods, firing rate and playback rate must be positive.")
    if values["readout_period_ms"] * values["firing_rate_hz"] > 1000 * (1 + 1e-12):
        raise ValueError(
            "Source firing rate exceeds detector capacity. Reduce firing rate or readout period; one frame is acquired per pulse."
        )
    if (
        isinstance(n_views, (bool, np.bool_))
        or not isinstance(n_views, (int, np.integer))
        or not 1 <= n_views <= MAX_RECIPE_VIEWS
    ):
        raise ValueError(f"Views must be between 1 and {MAX_RECIPE_VIEWS}.")
    times = np.arange(n_views, dtype=float) / values["firing_rate_hz"]
    return values, times


def _view_angles(parameters, count):
    acquisition = parameters.get("acquisition")
    if acquisition is not None:
        values, times = acquisition_timing(count, acquisition)
        return (
            np.deg2rad(parameters.get("orbit_start_deg", 0))
            + 2 * np.pi * times / values["rotation_period_s"]
        )
    span = parameters.get("scan_span_deg", 360)
    return np.deg2rad(parameters.get("orbit_start_deg", 0)) + np.linspace(
        0, np.deg2rad(span), count, endpoint=span < 360
    )


def _source_segments(waypoints, segment_views):
    """Validate connected waypoints or explicit (segment, endpoint, xyz) pairs."""
    points = np.asarray(waypoints, dtype=float)
    counts = np.asarray(segment_views)
    if points.ndim == 2 and points.shape[1] == 3 and len(points) >= 2:
        segments = np.stack((points[:-1], points[1:]), axis=1)
    elif points.ndim == 3 and points.shape[1:] == (2, 3) and len(points) >= 1:
        segments = points
    else:
        raise ValueError(
            "Source waypoints require at least two finite x, y, z triples or a list of [start, end] pairs."
        )
    if not np.isfinite(segments).all():
        raise ValueError("Source segment coordinates must be finite.")
    if counts.shape != (len(segments),) or counts.dtype.kind not in "iu" or np.any(counts < 2):
        raise ValueError("Enter one integer view count of at least 2 for each source segment.")
    shared = np.r_[False, np.all(segments[1:, 0] == segments[:-1, 1], axis=1)]
    total = sum(int(count) - int(skip) for count, skip in zip(counts, shared))
    if total > MAX_RECIPE_VIEWS:
        raise ValueError(f"Piecewise source exceeds the {MAX_RECIPE_VIEWS} view limit.")
    return segments, counts, shared


def _piecewise_source_positions(waypoints, segment_views):
    """Sample source segments without acquiring views across spatial gaps."""
    segments, counts, shared = _source_segments(waypoints, segment_views)
    return np.concatenate(
        [
            np.linspace(start, end, int(count))[int(skip) :]
            for (start, end), count, skip in zip(segments, counts, shared)
        ]
    )


def _with_source_beam(surface, specification, waypoints=None, segment_views=None):
    if specification is None:
        return surface
    from tigre.utilities.source_beam import SourceBeam
    from tigre.utilities.detector_geometry import DetectorSurfaceGeometry

    spec = dict(specification)
    mode = spec.pop("axis_mode", "Toward origin")
    roll = np.deg2rad(spec.pop("roll_deg", 0.0))
    if not np.isfinite(roll):
        raise ValueError("Beam roll must be finite")
    if "axis_directions" in spec:
        axis = np.asarray(spec.pop("axis_directions"), dtype=float)
        axis = np.broadcast_to(axis, surface.source_positions_mm.shape).copy()
    elif mode == "Toward origin":
        axis = -surface.source_positions_mm.astype(float)
    elif mode == "Segment normal":
        if waypoints is None:
            raise ValueError("Segment-normal beams require a linear or piecewise linear source")
        segments, counts, shared = _source_segments(waypoints, segment_views)
        axes = []
        for (start, end), count, skip in zip(segments, counts, shared):
            tangent = end - start
            if tangent[2] != 0 or np.linalg.norm(tangent[:2]) == 0:
                raise ValueError(
                    "Segment-normal beams require nonzero source segments in the xy plane"
                )
            normal = np.array([tangent[1], -tangent[0], 0.0])
            if normal @ (-(start + end) / 2) < 0:
                normal *= -1
            axes.extend([normal] * (int(count) - int(skip)))
        axis = np.asarray(axes)
    else:
        raise ValueError("Beam axis mode must be 'Toward origin' or 'Segment normal'")
    lengths = np.linalg.norm(axis, axis=1)
    if np.any(lengths == 0) or not np.isfinite(axis).all():
        raise ValueError("Cannot aim a beam with a zero or nonfinite axis")
    axis /= lengths[:, None]
    # Zero roll: wide axis lies in xy, narrow axis is approximately world z.
    reference = np.tile([0.0, 0.0, 1.0], (len(axis), 1))
    reference[np.abs(axis[:, 2]) > 0.99] = [0.0, 1.0, 0.0]
    wide = np.cross(reference, axis)
    wide /= np.linalg.norm(wide, axis=1)[:, None]
    wide = wide * np.cos(roll) + np.cross(axis, wide) * np.sin(roll)
    wide = spec.pop("wide_axis_directions", wide)
    beam = SourceBeam(axis, wide, **spec)
    return DetectorSurfaceGeometry(
        surface.source_positions_mm,
        surface.pixel_centers_mm,
        surface.valid_mask,
        module_slices=surface._module_slices,
        source_beam=beam,
    )


def build_modular_surface(
    layout,
    n_views,
    source_radius_mm,
    source_z_mm,
    module_shape,
    pixel_spacing_mm,
    detector_radius_mm,
    column_gap,
    n_modules,
    angular_span_deg,
    transverse_centers_mm,
    normal_offsets_mm,
    panel_angles_deg,
    z_offsets_mm,
    source_motion="Rotating source + detector",
    source_end_mm=(1000.0, 200.0, 0.0),
    source_start_mm=None,
    orbit_start_deg=0.0,
    scan_span_deg=360.0,
    panel_x_angles_deg=(0.0,),
    panel_y_angles_deg=(0.0,),
    cor_mm=0.0,
    acquisition=None,
    fixed_hardware=False,
    hardware_angles_deg=None,
    top_arm_rows=0,
    arm_clearance_mm=0.0,
    source_waypoints_mm=None,
    source_segment_views=None,
    source_beam=None,
):
    """Build a surface using distances in mm and angles in degrees.

    Panel x/y angles accept one value or one per panel, constant per scan.
    Rotation about each panel centre is Rz(gantry + panel_z) @ Ry(y) @ Rx(x).
    Cylindrical panel_z angles come from the arc layout. Rotating-mode COR
    translates the shared source and every panel by
    cor_mm * (-sin(gantry), cos(gantry), 0). It is one global scalar;
    use panel position offsets and tilts for individual panel alignment.
    Returned coordinates include these corrections; do not apply them again.
    Fixed panels face the initial linear-source position or the arc midpoint.
    source_waypoints_mm and source_segment_views define one or more linear
    segments on the same fixed detector, counting each shared junction once.
    Use a list of [start, end] pairs to leave unacquired gaps between segments.
    """
    if not isinstance(n_views, (int, np.integer)) or n_views < 1:
        raise ValueError("Number of views must be a positive integer.")
    if isinstance(n_views, (bool, np.bool_)):
        raise ValueError("Number of views must be a positive integer.")
    piecewise_source = None
    if source_waypoints_mm is not None or source_segment_views is not None:
        if (
            source_motion != "Fixed detector, moving source"
            or fixed_hardware
            or hardware_angles_deg is not None
        ):
            raise ValueError(
                "Piecewise source movement requires a fixed detector and moving source."
            )
        if acquisition is not None:
            raise ValueError("Piecewise source movement uses segment views; disable orbit timing.")
        piecewise_source = _piecewise_source_positions(source_waypoints_mm, source_segment_views)
        if len(piecewise_source) != n_views:
            raise ValueError("Views must equal the sum of segment views minus shared junctions.")
        source_start_mm, source_end_mm = piecewise_source[0], piecewise_source[-1]
    if fixed_hardware or hardware_angles_deg is not None:
        if acquisition is not None:
            raise ValueError("Object rotation uses explicit angles; disable orbit timing.")
        rows, columns = _checked_module_shape(module_shape)
        validate_surface_size(
            n_views,
            rows,
            columns,
            _panel_count_for_layout(layout, n_modules, transverse_centers_mm),
            column_gap,
        )
        import inspect

        arguments = locals()
        parameters = {
            key: arguments[key] for key in inspect.signature(build_modular_surface).parameters
        }
        parameters.update(
            n_views=1,
            fixed_hardware=False,
            hardware_angles_deg=None,
            source_beam=None,
            source_motion="Rotating source + detector",
        )
        initial = build_modular_surface(**parameters)
        from tigre.utilities.detector_geometry import DetectorSurfaceGeometry

        source = np.repeat(initial.source_positions_mm, n_views, axis=0)
        pixels = np.repeat(initial.pixel_centers_mm, n_views, axis=0)
        if not fixed_hardware:
            from types import SimpleNamespace
            from tigre.utilities.object_motion import object_angles, rotation_zyz

            angles = object_angles(SimpleNamespace(object_angles_deg=hardware_angles_deg), n_views)
            for view, pose in enumerate(angles):
                rotation = rotation_zyz(pose)
                source[view] = source[view] @ rotation.T
                pixels[view] = pixels[view] @ rotation.T
        result = DetectorSurfaceGeometry(
            source,
            pixels,
            np.repeat(initial.valid_mask, n_views, axis=0),
            module_slices=initial._module_slices,
        )
        if source_beam is not None and not fixed_hardware:
            # Define beam on the fixed hardware, then transform its directions
            # with the rig; do not re-aim it after object-frame transformation.
            beam_surface = _with_source_beam(
                DetectorSurfaceGeometry(
                    np.repeat(initial.source_positions_mm, n_views, axis=0),
                    np.repeat(initial.pixel_centers_mm, n_views, axis=0),
                ),
                source_beam,
            )
            rotations = np.array([rotation_zyz(pose) for pose in angles])
            return DetectorSurfaceGeometry(
                source,
                pixels,
                result.valid_mask,
                module_slices=result._module_slices,
                source_beam=beam_surface.source_beam.rotated(rotations),
            )
        return _with_source_beam(result, source_beam)
    rows, columns = _checked_module_shape(module_shape)
    if (
        isinstance(top_arm_rows, (bool, np.bool_))
        or not isinstance(top_arm_rows, (int, np.integer))
        or not 0 <= top_arm_rows < rows
    ):
        raise ValueError("Top arm rows must be an integer from 0 to module rows minus one.")
    panel_count = _panel_count_for_layout(layout, n_modules, transverse_centers_mm)
    validate_surface_size(n_views, rows, columns, panel_count, column_gap)
    if acquisition is not None:
        stationary_timing = (
            source_motion == "Fixed detector, moving source"
            and source_start_mm is not None
            and np.array_equal(np.asarray(source_start_mm), np.asarray(source_end_mm))
        )
        if source_motion != "Rotating source + detector" and not stationary_timing:
            raise ValueError(
                "Timing requires rotating hardware or an identical fixed-source start/end position."
            )
        acquisition_timing(n_views, acquisition)
    if source_motion in ("Rotating source + detector", "Fixed detector, arc source"):
        if not np.isfinite(source_radius_mm) or source_radius_mm <= 0:
            raise ValueError("Source radius must be a positive finite number.")
        if not np.isfinite(source_z_mm):
            raise ValueError("Source height must be finite.")
        if not np.isfinite(orbit_start_deg):
            raise ValueError("Orbit start angle must be finite.")
        if (
            not np.isfinite(scan_span_deg)
            or scan_span_deg <= 0
            or (acquisition is None and scan_span_deg > 360)
        ):
            raise ValueError(
                "Scan span must be positive and at most 360 degrees unless acquisition timing is enabled."
            )
        source_start = np.array([source_radius_mm, 0.0, source_z_mm], dtype=np.float64)
        angles = _view_angles(
            dict(
                orbit_start_deg=orbit_start_deg,
                scan_span_deg=scan_span_deg,
                acquisition=acquisition,
            ),
            int(n_views),
        )
        source = np.column_stack(
            (
                source_radius_mm * np.cos(angles),
                source_radius_mm * np.sin(angles),
                np.full(int(n_views), source_z_mm),
            )
        )
    elif source_motion == "Fixed detector, moving source":
        if source_start_mm is None:
            source_start_mm = (source_radius_mm, 0.0, source_z_mm)
        source_start = np.asarray(source_start_mm, dtype=np.float64)
        if source_start.shape != (3,) or not np.isfinite(source_start).all():
            raise ValueError("Source start position must contain finite x, y, z coordinates.")
        source_end = np.asarray(source_end_mm, dtype=np.float64)
        if source_end.shape != (3,) or not np.isfinite(source_end).all():
            raise ValueError("Source end position must contain finite x, y, z coordinates.")
        orientation_source = source_start.copy()
        if np.hypot(orientation_source[0], orientation_source[1]) == 0:
            orientation_source[0] = 1.0
        source = orientation_source[None, :]
    else:
        raise ValueError(
            "Source motion must be 'Rotating source + detector' or "
            "'Fixed detector, moving source'."
        )

    arc_source = None
    if source_motion == "Fixed detector, arc source":
        arc_source = source.copy()
        # Build the fixed panel facing the middle of the acquisition arc.
        middle = np.deg2rad(orbit_start_deg + scan_span_deg / 2)
        source = np.array(
            [[source_radius_mm * np.cos(middle), source_radius_mm * np.sin(middle), source_z_mm]]
        )

    if layout == "Cylindrical":
        detector = cylindrical_modular_detector_surface(
            source,
            module_shape=module_shape,
            pixel_spacing_mm=pixel_spacing_mm,
            detector_radius_mm=detector_radius_mm,
            n_modules=n_modules,
            angular_span_rad=np.deg2rad(angular_span_deg),
            column_gap=column_gap,
        )
    elif layout == "Piecewise linear":
        arrays = (
            np.asarray(normal_offsets_mm, dtype=np.float64),
            np.asarray(transverse_centers_mm, dtype=np.float64),
            np.asarray(z_offsets_mm, dtype=np.float64),
        )
        if not arrays[0].size or any(values.ndim != 1 for values in arrays):
            raise ValueError("Enter one normal, transverse, and z offset per panel.")
        if len({len(values) for values in arrays}) != 1:
            raise ValueError("Normal, transverse, and z offset lists must have equal lengths.")
        panel_angles = np.asarray(panel_angles_deg, dtype=np.float64)
        if panel_angles.ndim != 1 or len(panel_angles) != len(arrays[0]):
            raise ValueError("Enter one panel angle per panel.")
        centers = np.column_stack(arrays)
        detector = piecewise_linear_modular_detector_surface(
            source,
            module_shape=module_shape,
            pixel_spacing_mm=pixel_spacing_mm,
            detector_radius_mm=detector_radius_mm,
            module_centers_mm=centers,
            panel_angles_rad=np.deg2rad(panel_angles),
            column_gap=column_gap,
        )
    else:
        raise ValueError("Layout must be 'Cylindrical' or 'Piecewise linear'.")

    # Rotate each panel about its own centre, in its existing local frame.
    x_angles = _panel_tilts(panel_x_angles_deg, panel_count, "Panel x angles")
    y_angles = _panel_tilts(panel_y_angles_deg, panel_count, "Panel y angles")
    azimuths = np.arctan2(source[:, 1], source[:, 0])
    if layout == "Cylindrical":
        z_angles = (
            np.linspace(
                np.deg2rad(angular_span_deg) / 2, -np.deg2rad(angular_span_deg) / 2, panel_count
            )
            if panel_count > 1
            else np.zeros(1)
        )
    else:
        z_angles = np.deg2rad(panel_angles_deg)
    for index in range(panel_count):
        if x_angles[index] == 0 and y_angles[index] == 0:
            continue
        start = index * (columns + column_gap)
        grid = detector.pixel_centers_mm[:, :, start : start + columns]
        for view, azimuth in enumerate(azimuths):
            base = _panel_rotation(azimuth + z_angles[index], 0, 0)
            rotation = _panel_rotation(azimuth + z_angles[index], x_angles[index], y_angles[index])
            center = grid[view].mean(axis=(0, 1), dtype=np.float64)
            grid[view] = (grid[view] - center) @ base @ rotation.T + center

    if (
        not np.isfinite(arm_clearance_mm)
        or arm_clearance_mm < 0
        or (arm_clearance_mm and not top_arm_rows)
    ):
        raise ValueError("Arm clearance must be finite and nonnegative, with top arm rows enabled.")
    if top_arm_rows:
        lower = rows - top_arm_rows
        pitch = float(np.atleast_1d(pixel_spacing_mm)[0])
        facets = []
        for panel in range(panel_count):
            start = panel * (columns + column_gap)
            grid = detector.pixel_centers_mm[:, :, start : start + columns]
            for view, azimuth in enumerate(azimuths):
                rotation = _panel_rotation(
                    azimuth + z_angles[panel], x_angles[panel], y_angles[panel]
                )
                hinge = grid[view, lower - 1].astype(float) + rotation[:, 2] * pitch / 2
                distance = (np.arange(top_arm_rows) + 0.5) * pitch
                grid[view, lower:] = hinge[None, :, :] + distance[:, None, None] * rotation[:, 0]
                if arm_clearance_mm:
                    edge = hinge.mean(axis=0)
                    ray = edge - detector.source_positions_mm[view]
                    height = float(ray @ rotation[:, 2])
                    if height <= 0:
                        raise ValueError(
                            "Arm clearance requires the source below the joint in the panel frame."
                        )
                    grid[view, lower:] += ray * (arm_clearance_mm / height)

            facets.extend(
                [
                    (slice(0, lower), slice(start, start + columns)),
                    (slice(lower, rows), slice(start, start + columns)),
                ]
            )
        detector._module_slices = facets

    cor_mm = float(cor_mm)
    if not np.isfinite(cor_mm):
        raise ValueError("COR must be finite.")
    if source_motion != "Rotating source + detector" and cor_mm != 0:
        raise ValueError("COR is supported only for rotating source + detector.")
    if cor_mm:
        # Same translation as TIGRE's surface_from_geometry convention.
        shift = cor_mm * np.column_stack(
            (-np.sin(azimuths), np.cos(azimuths), np.zeros(len(source)))
        )
        detector.source_positions_mm += shift
        detector.pixel_centers_mm += shift[:, None, None, :]

    if source_motion == "Fixed detector, moving source":
        if piecewise_source is not None:
            from tigre.utilities.detector_geometry import DetectorSurfaceGeometry

            result = DetectorSurfaceGeometry(
                piecewise_source,
                np.broadcast_to(
                    detector.pixel_centers_mm, (n_views,) + detector.pixel_centers_mm.shape[1:]
                ).copy(),
                np.broadcast_to(
                    detector.valid_mask, (n_views,) + detector.valid_mask.shape[1:]
                ).copy(),
                module_slices=detector._module_slices,
            )
            return _with_source_beam(result, source_beam, source_waypoints_mm, source_segment_views)
        fixed = linear_source_surface(
            detector.pixel_centers_mm[0],
            source_start,
            source_end,
            int(n_views),
            valid_mask=detector.valid_mask[0],
        )
        fixed._module_slices = detector._module_slices
        return _with_source_beam(fixed, source_beam, [source_start, source_end], [int(n_views)])
    if arc_source is not None:
        from tigre.utilities.detector_geometry import DetectorSurfaceGeometry

        result = DetectorSurfaceGeometry(
            arc_source,
            np.broadcast_to(
                detector.pixel_centers_mm, (int(n_views),) + detector.pixel_centers_mm.shape[1:]
            ).copy(),
            np.broadcast_to(
                detector.valid_mask, (int(n_views),) + detector.valid_mask.shape[1:]
            ).copy(),
            module_slices=detector._module_slices,
        )
        return _with_source_beam(result, source_beam)
    return _with_source_beam(detector, source_beam)


def _panel_tilts(values, count, name):
    values = np.atleast_1d(np.asarray(values, dtype=np.float64))
    if values.ndim != 1 or values.size not in (1, count) or not np.isfinite(values).all():
        raise ValueError(f"{name} must contain one finite value or one per panel.")
    return np.deg2rad(np.broadcast_to(values, (count,)))


def _panel_rotation(z, x, y):
    """Local Rx, then Ry, then existing panel/gantry Rz (radians)."""
    cz, sz, cx, sx, cy, sy = np.cos(z), np.sin(z), np.cos(x), np.sin(x), np.cos(y), np.sin(y)
    return (
        np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
        @ np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        @ np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    )


def _checked_module_shape(module_shape):
    shape = np.asarray(module_shape)
    if (
        shape.shape != (2,)
        or shape.dtype.kind not in "iu"
        or shape.dtype.kind == "b"
        or np.any(shape <= 0)
    ):
        raise ValueError("Module rows and columns must be positive integers.")
    rows, columns = (int(value) for value in shape)
    if rows > MAX_MODULE_DIMENSION or columns > MAX_MODULE_DIMENSION:
        raise ValueError(f"Each module dimension is limited to {MAX_MODULE_DIMENSION} pixels.")
    return rows, columns


def _panel_count_for_layout(layout, n_modules, transverse_centers_mm):
    if layout == "Cylindrical":
        if (
            not isinstance(n_modules, (int, np.integer))
            or isinstance(n_modules, (bool, np.bool_))
            or n_modules < 1
        ):
            raise ValueError("Module count must be a positive integer.")
        count = int(n_modules)
    elif layout == "Piecewise linear":
        count = len(np.atleast_1d(transverse_centers_mm))
        if count < 1:
            raise ValueError("Enter at least one piecewise-linear panel.")
    else:
        raise ValueError("Layout must be 'Cylindrical' or 'Piecewise linear'.")
    if count > MAX_MODULES:
        raise ValueError(f"Recipe construction supports at most {MAX_MODULES} detector modules.")
    return count


def validate_surface_size(n_views, rows, columns, panel_count, column_gap):
    """Reject recipe inputs whose generated per-pixel geometry is too large."""
    if n_views > MAX_RECIPE_VIEWS:
        raise ValueError(f"Recipe construction supports at most {MAX_RECIPE_VIEWS} views.")
    if (
        not isinstance(column_gap, (int, np.integer))
        or isinstance(column_gap, (bool, np.bool_))
        or column_gap < 0
    ):
        raise ValueError("Gap columns must be a non-negative integer.")
    dense_columns = panel_count * columns + max(0, panel_count - 1) * column_gap
    total_pixels = n_views * rows * dense_columns
    if total_pixels > MAX_TOTAL_PIXELS:
        raise ValueError(
            f"This layout requires {total_pixels:,} detector samples; "
            f"the recipe limit is {MAX_TOTAL_PIXELS:,}. Reduce views or detector size."
        )


def check_tigre_compatibility(surface):
    """Check surface-projector validity and whether the geometry fits surface FDK."""
    source = np.asarray(surface.source_positions_mm)
    pixels = np.asarray(surface.pixel_centers_mm)
    mask = np.asarray(surface.valid_mask)
    if source.ndim != 2 or source.shape[1] != 3 or source.shape[0] < 1:
        raise ValueError("Source positions must have shape (views, 3).")
    if (
        pixels.ndim != 4
        or pixels.shape[0] != source.shape[0]
        or pixels.shape[-1] != 3
        or min(pixels.shape[1:3]) < 1
    ):
        raise ValueError("Detector pixels must have shape (views, rows, columns, 3).")
    if mask.shape != pixels.shape[:3] or mask.dtype.kind != "b":
        raise ValueError("Detector valid_mask must be boolean with shape (views, rows, columns).")
    if not np.isfinite(source).all() or not np.isfinite(pixels).all():
        raise ValueError("Source and detector coordinates must be finite.")
    if not mask.any():
        raise ValueError("The detector has no valid pixels.")

    report = {
        "iterative_supported": True,
        "iterative_message": (
            "Supported by Ax_surface/Atb_surface and surface iterative algorithms "
            "(cone-beam Siddon; GPU runtime/build not tested)."
        ),
        "fdk_supported": False,
        "fdk_message": "",
        "fdk_available": False,
    }
    try:
        from tigre.utilities._surface_geometry_checks import _analyse

        if source.shape[0] < 2:
            raise ValueError("FDK requires at least two views.")
        _analyse(surface)
    except (ValueError, NotImplementedError) as exc:
        report["fdk_message"] = (
            f"Not supported by the validated circular path: {exc} "
            "Surface FDK reconstruction is not included in this API contribution."
        )
    else:
        report["fdk_supported"] = True
        report["fdk_message"] = (
            "Geometry meets circular surface-FDK restrictions; reconstruction "
            "is not included in this API contribution. This is not an accuracy guarantee."
        )
    return report


# Compatibility aliases for pre-split scripts and the initial GUI release.
MAX_GUI_VIEWS = MAX_RECIPE_VIEWS
validate_gui_size = validate_surface_size
