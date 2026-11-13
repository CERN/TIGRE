"""Rigid object rotation about its centre, with ZYZ Euler orientations."""
from types import SimpleNamespace
import numpy as np


def object_angles(geometry, count):
    angles = np.atleast_1d(np.asarray(getattr(geometry, "object_angles_deg", [0.0]), dtype=float))
    if angles.ndim == 1:
        angles = np.column_stack((angles, np.zeros((len(angles), 2))))
    if (angles.ndim != 2 or angles.shape[1] != 3 or len(angles) not in (1, count)
            or not np.isfinite(angles).all()):
        raise ValueError("Object angles require one finite ZYZ triple or one triple per view.")
    return np.broadcast_to(angles, (count, 3))


def rotation_z(degrees):
    angle = np.deg2rad(degrees)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rotation_zyz(degrees):
    """Active object pose Rz(z1) @ Ry(y) @ Rz(z2), angles in degrees."""
    z1, y, z2 = degrees
    c, s = np.cos(np.deg2rad(y)), np.sin(np.deg2rad(y))
    return rotation_z(z1) @ np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]) @ rotation_z(z2)


def fixed_object_geometry(surface, geometry):
    """Map world coordinates to object coordinates; preserve projection ordering."""
    from .detector_geometry import DetectorSurfaceGeometry
    count = len(surface.source_positions_mm)
    angles = object_angles(geometry, count)
    centers = np.atleast_2d(geometry.offOrigin)[:, ::-1]
    if len(centers) not in (1, count):
        raise ValueError("Volume offsets must have one entry or one per view.")
    source = np.empty_like(surface.source_positions_mm)
    pixels = np.empty_like(surface.pixel_centers_mm)
    for view, angle in enumerate(angles):
        center = centers[0 if len(centers) == 1 else view]
        # Row-vector inverse rotation is multiplication by R, not R.T.
        rotation = rotation_zyz(angle)
        source[view] = (surface.source_positions_mm[view] - center) @ rotation
        pixels[view] = (surface.pixel_centers_mm[view] - center) @ rotation
    beam = getattr(surface, "source_beam", None)
    if beam is not None:
        beam = beam.rotated(np.array([rotation_zyz(a).T for a in angles]))
    transformed = DetectorSurfaceGeometry(source, pixels,
        getattr(surface, "detector_valid_mask", surface.valid_mask),
        module_slices=getattr(surface, "_module_slices", None), source_beam=beam)
    volume = SimpleNamespace(sVoxel=geometry.sVoxel, offOrigin=np.zeros(3))
    return transformed, volume
