"""Elliptical angular source contours and independent hard apertures.

Angles are FULL openings in degrees. Ellipses are defined in the tangent
plane normal to the beam axis. A flux contour is not a cutoff.
"""
import numpy as np


class SourceBeam:
    """One beam frame per view; angles/fractions may be shared or per-view.

    ``wide_axis_directions`` fixes ellipse roll in world coordinates. It must
    be perpendicular to the central axis. ``opening_angles_deg`` holds
    (wide, narrow) full openings at ``reference_fraction`` of central flux.
    An optional elliptical ``aperture_angles_deg`` clips rays independently.
    """

    def __init__(self, axis_directions, wide_axis_directions, opening_angles_deg,
                 reference_fraction=0.5, aperture_angles_deg=None):
        axis = np.atleast_2d(np.asarray(axis_directions, dtype=float))
        if axis.ndim != 2 or axis.shape[1] != 3 or not len(axis):
            raise ValueError("Beam axes must have shape (views, 3)")
        n = len(axis)
        def vectors(value):
            try:
                a = np.broadcast_to(np.asarray(value, dtype=float), (n, 3)).copy()
            except ValueError as exc:
                raise ValueError("Beam directions must contain one xyz triple per view") from exc
            lengths = np.linalg.norm(a, axis=1)
            if not np.isfinite(a).all() or np.any(lengths == 0):
                raise ValueError("Beam directions must be finite and nonzero")
            return a / lengths[:, None]
        self.axis_directions = vectors(axis)
        self.wide_axis_directions = vectors(wide_axis_directions)
        if np.any(np.abs(np.sum(self.axis_directions * self.wide_axis_directions, axis=1)) > 1e-6):
            raise ValueError("Beam wide axis must be perpendicular to central axis")
        def angles(value):
            try:
                a = np.broadcast_to(np.asarray(value, dtype=float), (n, 2)).copy()
            except ValueError as exc:
                raise ValueError("Beam angles require (wide, narrow) or (views, 2)") from exc
            if (not np.isfinite(a).all() or np.any(a <= 0) or np.any(a >= 180)
                    or np.any(a[:, 0] < a[:, 1])):
                raise ValueError("Full beam angles must satisfy 0 < narrow <= wide < 180 degrees")
            return a
        self.opening_angles_deg = angles(opening_angles_deg)
        self.reference_fraction = np.broadcast_to(np.asarray(reference_fraction, dtype=float), (n,)).copy()
        if (not np.isfinite(self.reference_fraction).all() or
                np.any((self.reference_fraction <= 0) | (self.reference_fraction >= 1))):
            raise ValueError("Beam reference fraction must be between 0 and 1")
        self.aperture_angles_deg = None if aperture_angles_deg is None else angles(aperture_angles_deg)

    @property
    def narrow_axis_directions(self):
        return np.cross(self.axis_directions, self.wide_axis_directions)

    def subset(self, indices):
        return SourceBeam(self.axis_directions[indices], self.wide_axis_directions[indices],
            self.opening_angles_deg[indices], self.reference_fraction[indices],
            None if self.aperture_angles_deg is None else self.aperture_angles_deg[indices])

    def rotated(self, rotations):
        """Return frames transformed by one active rotation per view."""
        return SourceBeam(np.einsum('vij,vj->vi', rotations, self.axis_directions),
            np.einsum('vij,vj->vi', rotations, self.wide_axis_directions),
            self.opening_angles_deg, self.reference_fraction, self.aperture_angles_deg)

    def _radius_squared(self, directions, view, angles):
        forward = directions @ self.axis_directions[view]
        wide = directions @ self.wide_axis_directions[view]
        narrow = directions @ self.narrow_axis_directions[view]
        tangents = np.tan(np.deg2rad(angles[view] / 2))
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            q = (wide / (forward * tangents[0])) ** 2 + (narrow / (forward * tangents[1])) ** 2
        return np.where(forward > 0, q, np.inf)

    def aperture_mask(self, sources, pixels):
        mask = np.ones(pixels.shape[:3], dtype=bool)
        if self.aperture_angles_deg is not None:
            for view in range(len(sources)):
                q = self._radius_squared(pixels[view].astype(float) - sources[view],
                                         view, self.aperture_angles_deg)
                mask[view] = q <= 1 + 1e-12
        return mask

    def relative_flux(self, directions, view=0):
        """Optional Gaussian-in-tangent-angle approximation, not measured data.

        The selected contour has exactly the reference fraction. This profile
        is NOT multiplied into log line integrals by Ax_surface/Atb_surface.
        Hard aperture is applied separately when requested.
        """
        directions = np.asarray(directions, dtype=float)
        q = self._radius_squared(directions, view, self.opening_angles_deg)
        flux = self.reference_fraction[view] ** q
        if self.aperture_angles_deg is not None:
            flux = np.where(self._radius_squared(directions, view, self.aperture_angles_deg) <= 1 + 1e-12,
                            flux, 0)
        return flux

    def outline(self, source, view, distance, aperture=False):
        """Sample a reference/aperture ellipse for beam-envelope visualization."""
        angles = self.aperture_angles_deg if aperture else self.opening_angles_deg
        if angles is None:
            return None
        theta = np.linspace(0, 2 * np.pi, 65)
        halfwidth = distance * np.tan(np.deg2rad(angles[view] / 2))
        return (source + distance * self.axis_directions[view]
                + halfwidth[0] * np.cos(theta)[:, None] * self.wide_axis_directions[view]
                + halfwidth[1] * np.sin(theta)[:, None] * self.narrow_axis_directions[view])
