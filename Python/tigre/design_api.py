"""Supported interface for external geometry-design applications.

Applications such as ``tigre-geometry-design-tool`` must import TIGRE geometry
building, file I/O, validation and rendering helpers from here rather than from
``tigre.utilities``.  Internal modules may be reorganised freely; this module
keeps the names below stable.

``API_VERSION`` is ``(major, minor)``.  The major number changes when a name is
removed or its behaviour is incompatible; the minor number changes when names
are added.  Clients should check it with :func:`require_api`.

Importing this module never imports tkinter or any desktop application.
"""

from .utilities.cargo_detector import cargo_module_layout
from .utilities.ct_physics import validate_transmission_geometry
from .utilities.detector_geometry import (
    DetectorSurfaceGeometry,
    linear_source_surface,
    surface_from_geometry,
)
from .utilities.detector_geometry_builder import (
    MAX_MODULE_DIMENSION,
    MAX_MODULES,
    MAX_RECIPE_VIEWS,
    MAX_TOTAL_PIXELS,
    _checked_module_shape as checked_module_shape,
    _panel_count_for_layout as panel_count_for_layout,
    _panel_rotation as panel_rotation,
    _panel_tilts as panel_tilts,
    _piecewise_source_positions as piecewise_source_positions,
    _source_segments as source_segments,
    _view_angles as view_angles,
    _with_source_beam as with_source_beam,
    acquisition_timing,
    build_modular_surface,
    check_tigre_compatibility,
    validate_surface_size,
)
from .utilities.detector_geometry_remap import remap_projections_for_fdk
from .utilities.geometry_io import (
    _volume as volume_from_mapping,
    load_geometry_archive,
    load_geometry_file,
    save_geometry_yaml,
    surface_archive_arrays,
)
from .utilities.object_motion import fixed_object_geometry, object_angles, rotation_zyz
from .utilities.sample_loader import load_head_phantom
from .utilities.visualization.geometry_markers import draw_world_axes
from .utilities.visualization.surface_renderer import (
    SurfaceRenderer,
    _volume_box as volume_box,
)

API_VERSION = (1, 1)

__all__ = [
    "API_VERSION",
    "require_api",
    "DetectorSurfaceGeometry",
    "MAX_MODULES",
    "MAX_MODULE_DIMENSION",
    "MAX_RECIPE_VIEWS",
    "MAX_TOTAL_PIXELS",
    "SurfaceRenderer",
    "acquisition_timing",
    "build_modular_surface",
    "cargo_module_layout",
    "check_tigre_compatibility",
    "checked_module_shape",
    "draw_world_axes",
    "fixed_object_geometry",
    "linear_source_surface",
    "load_geometry_archive",
    "load_geometry_file",
    "load_head_phantom",
    "object_angles",
    "panel_count_for_layout",
    "panel_rotation",
    "panel_tilts",
    "piecewise_source_positions",
    "remap_projections_for_fdk",
    "rotation_zyz",
    "save_geometry_yaml",
    "source_segments",
    "surface_archive_arrays",
    "surface_from_geometry",
    "validate_surface_size",
    "validate_transmission_geometry",
    "view_angles",
    "volume_box",
    "volume_from_mapping",
    "with_source_beam",
]


def require_api(major, minimum_minor=0):
    """Raise ``RuntimeError`` unless this API is ``major`` with at least ``minimum_minor``."""
    found_major, found_minor = API_VERSION
    if found_major != major or found_minor < minimum_minor:
        raise RuntimeError(
            "TIGRE design API {}.{} is installed, but {}.{} or newer (same major) is required.".format(
                found_major, found_minor, major, minimum_minor
            )
        )


# Additive application services; all 1.0 exports above remain available.
from .application_api import (AlgorithmInfo, DeviceInfo, ReconstructionProgress, ReconstructionResult, ReconstructionCancelled, BackendExecutionError, available_algorithms, backend_capabilities, available_devices, create_operator, forward_project, reconstruct)
from .application_api import __all__ as _application_names
__all__ += _application_names
del _application_names
