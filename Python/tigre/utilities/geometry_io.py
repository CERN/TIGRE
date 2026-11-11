"""Versioned YAML settings and exact numeric surface archives.

YAML uses mm and degrees for modular builder settings. Explicit pixel arrays
remain in NPZ companions, loaded without pickle. No Tk window is required.
"""
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Optional
import inspect
import uuid
import tempfile
import numpy as np


def surface_archive_arrays(surface, geometry=None):
    """Numeric archive payload including optional module and volume metadata."""
    arrays = {name: getattr(surface, name) for name in (
        "source_positions_mm", "pixel_centers_mm", "valid_mask")}
    modules = getattr(surface, "_module_slices", None)
    beam = getattr(surface, "source_beam", None)
    if beam is not None:
        arrays["detector_valid_mask"] = surface.detector_valid_mask
        for name in ("axis_directions", "wide_axis_directions", "opening_angles_deg", "reference_fraction", "aperture_angles_deg"):
            value = getattr(beam, name)
            if value is not None:
                arrays["beam_" + name] = value
    if modules:
        arrays["module_slices"] = np.array([
            [r.start, r.stop, c.start, c.stop] for r, c in modules], dtype=np.int64)
    if geometry is not None:
        from tigre.utilities.visualization.surface_renderer import _volume_box
        origin, size = _volume_box(geometry)
        arrays["sVoxel"] = size[::-1]
        arrays["offOrigin"] = np.asarray(geometry.offOrigin, dtype=float)
        _validate_volume_views(geometry, len(surface.source_positions_mm))
        if hasattr(geometry, "object_angles_deg"):
            arrays["object_angles_deg"] = np.asarray(geometry.object_angles_deg)
        if hasattr(geometry, "simulation_mode"):
            arrays["simulation_mode"] = geometry.simulation_mode
        if hasattr(geometry, "motion_mode"):
            arrays["volume_motion_mode"] = geometry.motion_mode
        if hasattr(geometry, "nVoxel"):
            arrays["nVoxel"] = _volume(dict(sVoxel=size[::-1],
                nVoxel=geometry.nVoxel)).nVoxel
        if hasattr(geometry, "acquisition"):
            from .detector_geometry_builder import acquisition_timing
            timing, times = acquisition_timing(len(surface.source_positions_mm), geometry.acquisition)
            arrays["view_times_s"] = times
            arrays.update({"acquisition_" + name: value for name, value in timing.items()})
    return arrays


def load_geometry_archive(path):
    """Load explicit surfaces or standard TIGRE parameters from numeric NPZ."""
    from types import SimpleNamespace
    from tigre.utilities.detector_geometry import DetectorSurfaceGeometry, surface_from_geometry
    with np.load(path, allow_pickle=False) as archive:
        keys = {"source_positions_mm", "pixel_centers_mm", "valid_mask", "module_slices",
            "sVoxel", "offOrigin", "DSO", "DSD", "nDetector", "dDetector", "angles",
            "offDetector", "rotDetector", "COR", "nVoxel", "mode"}
        timing_names = ("readout_period_ms", "firing_rate_hz", "rotation_period_s",
                        "table_speed_mm_s", "playback_rate")
        keys.update("acquisition_" + name for name in timing_names)
        keys.update(("view_times_s", "volume_motion_mode", "object_angles_deg", "simulation_mode"))
        beam_names = ("axis_directions", "wide_axis_directions", "opening_angles_deg", "reference_fraction", "aperture_angles_deg")
        keys.update("beam_" + name for name in beam_names)
        keys.add("detector_valid_mask")
        arrays = {key: archive[key] for key in archive.files if key in keys}
    if "source_positions_mm" in arrays:
        modules = arrays.get("module_slices")
        if modules is not None:
            if modules.ndim != 2 or modules.shape[1] != 4 or modules.dtype.kind not in "iu":
                raise ValueError("module_slices must be an integer (modules, 4) array.")
            modules = [(slice(int(a), int(b)), slice(int(c), int(d))) for a, b, c, d in modules]
        beam_values = {name: arrays["beam_" + name] for name in beam_names if "beam_" + name in arrays}
        beam = None
        if beam_values:
            from .source_beam import SourceBeam
            if not {"axis_directions", "wide_axis_directions", "opening_angles_deg", "reference_fraction"} <= beam_values.keys():
                raise ValueError("Incomplete source beam metadata in geometry archive")
            beam = SourceBeam(**beam_values)
        surface = DetectorSurfaceGeometry(arrays["source_positions_mm"],
            arrays["pixel_centers_mm"], arrays.get("detector_valid_mask", arrays.get("valid_mask")),
            module_slices=modules, source_beam=beam)
        if beam is not None and "valid_mask" in arrays and not np.array_equal(surface.valid_mask, arrays["valid_mask"]):
            raise ValueError("Saved effective mask does not match beam aperture and detector mask")
        geometry = None
        if "sVoxel" in arrays:
            size = np.asarray(arrays["sVoxel"], dtype=float)
            offset = np.asarray(arrays.get("offOrigin", np.zeros(3)), dtype=float)
            if size.shape != (3,) or np.any(size <= 0) or not np.isfinite(size).all():
                raise ValueError("sVoxel must contain three positive finite values.")
            if offset.ndim not in (1, 2) or offset.shape[-1] != 3 or not np.isfinite(offset).all():
                raise ValueError("offOrigin must contain three finite values.")
            geometry = _volume(dict(sVoxel=size, offOrigin=offset,
                **({"nVoxel": arrays["nVoxel"]} if "nVoxel" in arrays else {}),
                **({"simulation_mode": arrays["simulation_mode"].item()} if "simulation_mode" in arrays else {}),
                **({"object_angles_deg": arrays["object_angles_deg"]} if "object_angles_deg" in arrays else {}),
                **({"motion_mode": arrays["volume_motion_mode"].item()} if "volume_motion_mode" in arrays else {})))
            _validate_volume_views(geometry, len(surface.source_positions_mm))
            if any("acquisition_" + name in arrays for name in timing_names):
                from .detector_geometry_builder import acquisition_timing
                timing = {name: arrays["acquisition_" + name].item() for name in timing_names
                          if "acquisition_" + name in arrays}
                geometry.acquisition, geometry.view_times_s = acquisition_timing(
                    len(surface.source_positions_mm), timing)
                if "view_times_s" in arrays and not np.array_equal(arrays["view_times_s"], geometry.view_times_s):
                    raise ValueError("Saved view timestamps do not match acquisition timing.")
        return geometry, surface
    required = ("DSO", "DSD", "nDetector", "dDetector", "angles")
    if any(name not in arrays for name in required):
        raise ValueError("Archive needs source/pixel arrays, or DSO, DSD, nDetector, dDetector and angles.")
    if "mode" in arrays and arrays["mode"].item() != "cone":
        raise ValueError("Standard geometry inspection supports cone mode; supply an explicit surface otherwise.")
    import tigre
    geometry = tigre.geometry_default(high_resolution=False)
    for name in ("DSO", "DSD", "nDetector", "dDetector", "offDetector", "rotDetector", "COR",
                 "nVoxel", "sVoxel", "offOrigin"):
        if name in arrays:
            value = arrays[name]
            setattr(geometry, name, value.item() if value.ndim == 0 else value)
    geometry.sDetector = geometry.nDetector * geometry.dDetector
    geometry.dVoxel = geometry.sVoxel / geometry.nVoxel
    return geometry, surface_from_geometry(geometry, arrays["angles"])



@dataclass
class GeometryDocument:
    geometry: object
    surface: object
    settings: Optional[dict] = None


def _yaml():
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("YAML geometry files require PyYAML: install it in your TIGRE environment.") from exc
    return yaml


def _plain(value):
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    return value


def _settings(values):
    from .detector_geometry_builder import build_modular_surface
    if not isinstance(values, dict):
        raise ValueError("builder must be a mapping of modular geometry settings.")
    try:
        bound = inspect.signature(build_modular_surface).bind(**values)
    except TypeError as exc:
        raise ValueError(f"Invalid builder settings: {exc}") from exc
    bound.apply_defaults()
    return _plain(dict(bound.arguments))


def _volume(values):
    if values is None:
        return None
    if not isinstance(values, dict) or set(values) - {"sVoxel", "offOrigin", "nVoxel", "motion_mode", "object_angles_deg", "simulation_mode"} or "sVoxel" not in values:
        raise ValueError("volume requires sVoxel and optional offOrigin, in TIGRE z/y/x order.")
    geo = SimpleNamespace(sVoxel=np.asarray(values["sVoxel"], dtype=float),
                          offOrigin=np.asarray(values.get("offOrigin", [0, 0, 0]), dtype=float))
    from .visualization.surface_renderer import _volume_box
    _volume_box(geo)
    if "object_angles_deg" in values:
        geo.object_angles_deg = np.atleast_1d(np.asarray(values["object_angles_deg"], dtype=float))
        from .object_motion import object_angles
        object_angles(geo, len(geo.object_angles_deg))
    if "simulation_mode" in values:
        if values["simulation_mode"] not in ("ct", "dr"):
            raise ValueError("simulation_mode must be ct or dr.")
        geo.simulation_mode = values["simulation_mode"]
    if "motion_mode" in values:
        if values["motion_mode"] not in ("static", "speed", "endpoints"):
            raise ValueError("volume motion_mode must be static, speed or endpoints.")
        geo.motion_mode = values["motion_mode"]
    if "nVoxel" in values:
        counts = np.asarray(values["nVoxel"], dtype=float)
        if (counts.shape != (3,) or not np.isfinite(counts).all()
                or np.any(counts <= 0) or np.any(counts != np.floor(counts))
                or np.any(counts >= np.iinfo(np.int64).max)):
            raise ValueError("nVoxel must contain three positive integers in z/y/x order.")
        geo.nVoxel = counts.astype(np.int64)
        geo.dVoxel = geo.sVoxel / geo.nVoxel
    return geo


def _validate_volume_views(geo, n_views):
    if geo is not None:
        from .object_motion import object_angles
        object_angles(geo, n_views)
    if geo is not None and len(np.atleast_2d(geo.offOrigin)) not in (1, n_views):
        raise ValueError("offOrigin must contain one offset or one per view.")


def _write_yaml(path, text):
    """Replace a recipe only after the complete new file has been written."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=path.name + ".", suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        temporary.replace(path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def load_geometry_file(path):
    """Return geometry, exact/generated surface, and editable YAML settings.

    Version 1 YAML is either ``kind: modular`` with ``builder`` settings or
    ``kind: surface`` with a relative/local ``surface_file`` NPZ reference.
    """
    path = Path(path)
    if path.suffix.lower() == ".npz":
        geo, surface = load_geometry_archive(path)
        return GeometryDocument(geo, surface)
    if path.suffix.lower() not in (".yaml", ".yml"):
        raise ValueError("Geometry files must be .yaml, .yml or .npz.")
    yaml = _yaml()
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid geometry YAML: {exc}") from exc
    if not isinstance(data, dict) or data.get("schema") != "tigre.geometry":
        raise ValueError("Expected a tigre.geometry YAML document.")
    if type(data.get("version")) is not int or data["version"] != 1:
        raise ValueError("Unsupported geometry YAML version; expected 1.")
    if data.get("units") != {"length": "mm", "angle": "degree"}:
        raise ValueError("Geometry YAML units must be length: mm and angle: degree.")
    if data.get("rotation_convention") != "Rz(gantry+z) @ Ry(y) @ Rx(x)":
        raise ValueError("Unsupported rotation convention.")
    kind = data.get("kind")
    allowed = {"schema", "version", "units", "rotation_convention", "kind", "volume"}
    allowed.add("builder" if kind == "modular" else "surface_file")
    if set(data) - allowed:
        raise ValueError(f"Unknown YAML fields: {sorted(set(data) - allowed)}")
    geo = _volume(data.get("volume"))
    if kind == "modular":
        from .detector_geometry_builder import build_modular_surface
        settings = _settings(data.get("builder"))
        surface = build_modular_surface(**settings)
        _validate_volume_views(geo, len(surface.source_positions_mm))
        if geo is not None and settings.get("acquisition") is not None:
            from .detector_geometry_builder import acquisition_timing
            geo.acquisition, geo.view_times_s = acquisition_timing(
                settings["n_views"], settings["acquisition"])
        return GeometryDocument(geo, surface, settings)
    if kind == "surface":
        reference = data.get("surface_file")
        if not isinstance(reference, str) or Path(reference).suffix.lower() != ".npz":
            raise ValueError("surface_file must reference a numeric .npz archive.")
        stored_geo, surface = load_geometry_archive(path.parent / reference)
        _validate_volume_views(geo, len(surface.source_positions_mm))
        if geo is not None and hasattr(stored_geo, "acquisition"):
            geo.acquisition = stored_geo.acquisition
            geo.view_times_s = stored_geo.view_times_s
        return GeometryDocument(geo if geo is not None else stored_geo, surface)
    raise ValueError("Geometry YAML kind must be modular or surface.")


def save_geometry_yaml(path, *, settings=None, surface=None, geometry=None):
    """Save editable settings or an exact surface with a unique NPZ companion.

    Builder settings are validated by constructing their surface. Inspection
    surfaces are never fitted to an approximate editable geometry.
    """
    path = Path(path)
    if path.suffix.lower() not in (".yaml", ".yml"):
        raise ValueError("Geometry settings filename must end in .yaml or .yml.")
    yaml = _yaml()
    data = dict(schema="tigre.geometry", version=1,
                units={"length": "mm", "angle": "degree"},
                rotation_convention="Rz(gantry+z) @ Ry(y) @ Rx(x)")
    if geometry is not None:
        from .visualization.surface_renderer import _volume_box
        origin, size = _volume_box(geometry)
        data["volume"] = dict(sVoxel=size[::-1].tolist(),
                              offOrigin=np.asarray(geometry.offOrigin).tolist())
        if hasattr(geometry, "simulation_mode"):
            data["volume"]["simulation_mode"] = geometry.simulation_mode
        if hasattr(geometry, "object_angles_deg"):
            data["volume"]["object_angles_deg"] = np.asarray(geometry.object_angles_deg).tolist()
        if hasattr(geometry, "nVoxel"):
            data["volume"]["nVoxel"] = _volume(dict(
                **data["volume"], nVoxel=geometry.nVoxel)).nVoxel.tolist()
        if hasattr(geometry, "motion_mode"):
            data["volume"]["motion_mode"] = _volume(dict(
                **data["volume"], motion_mode=geometry.motion_mode)).motion_mode
    if settings is not None:
        from .detector_geometry_builder import build_modular_surface
        values = _settings(settings)
        built = build_modular_surface(**values)
        _validate_volume_views(geometry, len(built.source_positions_mm))
        data.update(kind="modular", builder=values)
    elif surface is not None:
        companion = path.with_name(path.stem + ".surface-" + uuid.uuid4().hex[:8] + ".npz")
        data.update(kind="surface", surface_file=companion.name)
        payload = surface_archive_arrays(surface, geometry)
        text = yaml.safe_dump(data, sort_keys=False)
        try:
            np.savez_compressed(companion, **payload)
            _write_yaml(path, text)
        except BaseException:
            if companion.exists():
                companion.unlink()
            raise
        return
    else:
        raise ValueError("Supply builder settings or an explicit detector surface.")
    _write_yaml(path, yaml.safe_dump(data, sort_keys=False))
