"""Shared detector-surface scene for script animation and GUI inspection."""
import numpy as np
import mpl_toolkits.mplot3d.art3d as art3d

def _detector_runs(surface):
    """Modules (or, without module info, contiguous valid column runs) as (start, stop)."""
    cols = surface.valid_mask.any(axis=(0, 1))
    slices = getattr(surface, "_module_slices", None)
    if slices:
        return [(c.start, c.stop) for _, c in slices if cols[c.start:c.stop].any()]
    runs, start = [], None
    for c, ok in enumerate(np.append(cols, False)):
        if ok and start is None:
            start = c
        elif not ok and start is not None:
            runs.append((start, c))
            start = None
    return runs


def _surface_mesh(surface, view, max_rows=24, max_columns=48):
    """Sample cells within module bounds; never span a masked pixel or module join.

    Vertices are measured pixel centres, not inferred physical pixel edges.
    Single-row/column detectors are represented by line segments.
    """
    pixels = surface.pixel_centers_mm[view]
    mask = surface.valid_mask[view]
    modules = getattr(surface, "_module_slices", None)
    if not modules:
        modules = [(slice(0, mask.shape[0]), slice(c0, c1))
                   for c0, c1 in _detector_runs(surface)]
    faces, segments = [], []
    for rows, columns in modules:
        r0, r1, c0, c1 = rows.start, rows.stop, columns.start, columns.stop
        rr = np.unique(np.linspace(r0, r1 - 1, min(r1 - r0, max_rows), dtype=int))
        cc = np.unique(np.linspace(c0, c1 - 1, min(c1 - c0, max_columns), dtype=int))
        if len(rr) == 1 or len(cc) == 1:
            indices = ([(rr[0], c) for c in cc] if len(rr) == 1
                       else [(r, cc[0]) for r in rr])
            for (ra, ca), (rb, cb) in zip(indices, indices[1:]):
                if mask[ra:rb + 1, ca:cb + 1].all():
                    segments.append([pixels[ra, ca], pixels[rb, cb]])
        else:
            for ra, rb in zip(rr, rr[1:]):
                for ca, cb in zip(cc, cc[1:]):
                    if mask[ra:rb + 1, ca:cb + 1].all():
                        # Triangles also represent curved/non-planar grids reliably.
                        p00, p01 = pixels[ra, ca], pixels[ra, cb]
                        p10, p11 = pixels[rb, ca], pixels[rb, cb]
                        faces.extend([[p00, p01, p11], [p00, p11, p10]])
    return faces, segments


def _volume_box(geo):
    """Validate volume metadata and return initial world xyz bounds."""
    size = np.asarray(geo.sVoxel, dtype=float)
    offsets = np.atleast_2d(np.asarray(geo.offOrigin, dtype=float))
    if size.shape != (3,) or not np.isfinite(size).all() or np.any(size <= 0):
        raise ValueError("sVoxel must contain three positive finite values.")
    if offsets.ndim != 2 or offsets.shape[1] != 3 or not offsets.size or not np.isfinite(offsets).all():
        raise ValueError("offOrigin must contain finite three-coordinate offsets.")
    return offsets[0][::-1], size[::-1]


class SurfaceRenderer:
    """Build a scene on existing 3D axes and update its artists by view index.

    ``geo`` optionally supplies the volume. ``panel_vertices`` supplies physical
    panel outlines for builder geometry; arbitrary loaded surfaces use a mesh.
    Coordinates are displayed as supplied, with no extra rotation/COR applied.
    """

    def __init__(self, ax, surface, geo=None, pos=0, max_points=4000,
                 world_axes=True, panel_vertices=None):
        src = np.asarray(surface.source_positions_mm)
        pix = np.asarray(surface.pixel_centers_mm)
        mask = np.asarray(surface.valid_mask)
        n_views = len(src)
        pos = 0 if pos is None or pos < 0 or pos >= n_views else int(pos)
        if not mask.any():
            raise ValueError("the detector surface has no valid pixels")

        origin, size = (np.zeros(3), np.zeros(3)) if geo is None else _volume_box(geo)
        from ..object_motion import object_angles, rotation_zyz
        angles = object_angles(geo, n_views)
        volume_mesh = None
        origins = np.zeros((1, 3)) if geo is None else np.atleast_2d(geo.offOrigin)[:, ::-1]
        if len(origins) not in (1, n_views):
            raise ValueError("offOrigin must contain one offset or one offset per view.")
        if geo is not None:

            volume_mesh = art3d.Poly3DCollection(calCube(origin, size), color="c", alpha=0.2)
            ax.add_collection3d(volume_mesh)
        orientation = None
        if hasattr(geo, "object_angles_deg"):
            orientation, = ax.plot([], [], [], color="teal", linewidth=3, label="Object +x")
        if world_axes:
            from .geometry_markers import draw_world_axes
            length = (max(float(np.min(size)) / 2, 1) if geo is not None
                      else max(float(np.linalg.norm(src, axis=1).max()), 1) * 0.12)
            draw_world_axes(ax, length)

        ax.plot3D(*src.T, color="crimson", linestyle="--", alpha=0.5, label="Source path")
        source = ax.scatter(*src[pos], color="crimson", s=35, label="Source")
        source_label = ax.text(*src[pos], "S")
        mesh = art3d.Poly3DCollection([], facecolors="cornflowerblue", edgecolors="navy",
                                     linewidths=0.2, alpha=0.35)
        lines = art3d.Line3DCollection([], colors="navy", linewidths=1.5)
        ax.add_collection3d(mesh, autolim=False)
        ax.add_collection3d(lines, autolim=False)
        cloud = ax.scatter([], [], [], s=1, color="navy", alpha=0.6)
        central_ray, = ax.plot([], [], [], color="pink")
        beam = getattr(surface, "source_beam", None)
        beam_lines = []
        if beam is not None:
            for label, color in (("Beam flux contour", "darkorange"), ("Hard aperture", "purple")):
                if label == "Hard aperture" and beam.aperture_angles_deg is None:
                    continue
                beam_lines.append(ax.plot([], [], [], color=color, linewidth=1, alpha=0.2, label=label)[0])
        beam_distances = np.empty(n_views) if beam is not None else None
        title = ax.set_title("")

        # Scan frame-wise to avoid copying a large acquisition into float64.
        lo = np.minimum(src.min(axis=0), np.minimum(origin - size / 2, 0))
        hi = np.maximum(src.max(axis=0), np.maximum(origin + size / 2, 0))
        extent = size / 2
        if np.any(angles):
            extent = np.full(3, np.linalg.norm(size) / 2)
        lo = np.minimum(lo, origins.min(axis=0) - extent)
        hi = np.maximum(hi, origins.max(axis=0) + extent)
        module_count = len(_detector_runs(surface))
        for frame in range(n_views):
            valid = pix[frame][mask[frame]]
            if len(valid):
                lo, hi = np.minimum(lo, valid.min(axis=0)), np.maximum(hi, valid.max(axis=0))
            if beam is not None:
                # Include physical detector pixels even outside the aperture.
                detector_mask = getattr(surface, "detector_valid_mask", mask)[frame]
                physical = pix[frame][detector_mask]
                if len(physical):
                    lo = np.minimum(lo, physical.min(axis=0))
                    hi = np.maximum(hi, physical.max(axis=0))
                axis = beam.axis_directions[frame]
                axial = pix[frame] @ axis - src[frame] @ axis
                distance = max(float(np.max(axial, where=detector_mask, initial=0)) * 1.05, 1.)
                beam_distances[frame] = distance
                # Beam outlines may extend beyond the physical scene limits;
                # enabling the beam must not zoom out to fit its profile.
        span = np.maximum(hi - lo, 1)
        # Leave vertical breathing room for thin detectors, preserving mm scale.
        middle_z = (lo[2] + hi[2]) / 2
        span[2] = max(span[2], 0.5 * max(span[0], span[1]))
        lo[2], hi[2] = middle_z - span[2] / 2, middle_z + span[2] / 2
        padding = span * 0.04
        for setter, low, high in zip((ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d), lo-padding, hi+padding):
            setter(low, high)
        ax.set_box_aspect(span)
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_zlabel("Z (mm)")

        def update(frame):
            if volume_mesh is not None:
                center = origins[0 if len(origins) == 1 else frame]
                faces = np.asarray(calCube(np.zeros(3), size))
                volume_mesh.set_verts(faces @ rotation_zyz(angles[frame]).T + center)
                if orientation is not None:
                    tip = center + rotation_zyz(angles[frame]) @ np.array([size[0] / 2, 0, 0])
                    orientation.set_data_3d(*np.stack((center, tip)).T)
            source._offsets3d = tuple(src[frame, :, None])
            source_label.set_position_3d(src[frame])
            if beam is not None:
                distance = beam_distances[frame]
                for aperture, line in enumerate(beam_lines):
                    ring = beam.outline(src[frame], frame, distance, aperture=bool(aperture))
                    if ring is None:
                        line.set_data_3d([], [], [])
                    else:
                        # Closed ellipse plus four spokes show the cone orientation.
                        points = list(ring)
                        for index in (0, 16, 32, 48):
                            points.extend(([np.nan] * 3, src[frame], ring[index]))
                        line.set_data_3d(*np.asarray(points).T)
            if panel_vertices is None:
                faces, segments = _surface_mesh(surface, frame)
            else:
                faces, segments = panel_vertices(frame), []
            mesh.set_verts(faces)
            lines.set_segments(segments)
            valid = pix[frame][mask[frame]]
            stride = max(1, int(np.ceil(len(valid) / max_points)))
            samples = valid[::stride]
            cloud._offsets3d = tuple(samples.T)
            if len(valid):
                centre = valid.mean(axis=0, dtype=np.float64)
                central_ray.set_data_3d(*np.stack((src[frame], centre)).T)
            else:
                central_ray.set_data_3d([], [], [])
            title.set_text(f"Detector surface geometry - view {frame}/{n_views}\n"
                           f"{module_count} detector module(s), {len(valid)} valid pixels")
            artists = (source, source_label, mesh, lines, cloud, central_ray, title)
            return artists + tuple(beam_lines) + ((volume_mesh,) if volume_mesh is not None else ())

        self.surface = surface
        self.mesh = mesh
        self.source = source
        self.volume = volume_mesh
        self.update = update
        self.update(pos)


def calCube(center, size):
    """Six volume-box faces without importing an animation or UI module."""
    corners=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                      [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],dtype=float)
    vertices=corners*np.asarray(size)/2+np.asarray(center)
    return vertices[np.array([[1,2,6,5],[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[0,3,7,4]])]
