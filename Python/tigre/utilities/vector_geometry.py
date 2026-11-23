"""Convert between TIGRE cone-beam geometries and per-view geometry vectors.

The vector layout is the one used by ASTRA's ``cone_vec`` geometry: one row of
12 numbers per view,

    ``(srcX, srcY, srcZ, dX, dY, dZ, uX, uY, uZ, vX, vY, vZ)``

where ``src`` is the source position, ``d`` the centre of the detector, ``u``
the vector from detector pixel ``(row 0, column 0)`` to ``(row 0, column 1)``
and ``v`` the vector from ``(row 0, column 0)`` to ``(row 1, column 0)``.
Row and column refer to TIGRE projection data, which has shape
``(n_views, nDetector[0], nDetector[1])``, i.e. ``(views, rows, columns)``.

All vectors are in TIGRE's world frame, in the same length unit as ``geo``
(normally mm): the image is fixed, its centre is at the origin, and its
voxel indices increase along +x, +y and +z. Only the row layout is borrowed
from ASTRA. Data and volume array conventions differ between the toolboxes,
so check them before reusing ASTRA vectors directly.

TIGRE describes each view with ``DSO``, ``DSD``, three ZYZ Euler ``angles``,
``offDetector``, ``rotDetector``, ``offOrigin`` and ``COR``. This is
redundant: every vector geometry (any rigid placement of source and detector
relative to the image) has many TIGRE descriptions. ``from_vector_to_TIGRE``
returns one of them, with ``offOrigin`` and ``COR`` set to zero. Vector
geometries that no TIGRE geometry can describe (detector axes that are not
perpendicular, or a pixel size that changes between views) are refused
rather than approximated.
"""

import copy

import numpy as np

__all__ = ["from_vector_to_TIGRE", "from_TIGRE_to_vector"]

_REFERENCES = ("center", "first_pixel")


def from_TIGRE_to_vector(geo, angles=None, detector_reference="center"):
    """Per-view source and detector vectors of a TIGRE cone-beam geometry.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry. ``DSO``, ``DSD``, ``offOrigin``,
        ``offDetector``, ``rotDetector`` and ``COR`` may be given once or per
        view. It is not modified.
    angles : numpy.ndarray, optional
        Projection angles, shape ``(N,)`` or ``(N, 3)``, in radians. Default
        ``geo.angles``.
    detector_reference : {"center", "first_pixel"}, optional
        Whether columns 3-5 hold the centre of the detector (default, as in
        ASTRA) or the centre of pixel ``(row 0, column 0)``.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(N, 12)``, ``float64``. See the module documentation
        for the layout.

    See Also
    --------
    from_vector_to_TIGRE
    """
    _check_reference(detector_reference)
    if getattr(geo, "mode", "cone") not in (None, "cone"):
        raise ValueError(
            "from_TIGRE_to_vector supports cone-beam geometries only, got mode=%r" % geo.mode
        )
    if angles is None:
        angles = getattr(geo, "angles", None)
        if angles is None:
            raise ValueError("No angles given and geo.angles is not set")
    angles = np.asarray(angles, dtype=np.float64)
    if angles.ndim == 1:
        angles = np.stack([angles, np.zeros_like(angles), np.zeros_like(angles)], axis=1)
    if angles.ndim != 2 or angles.shape[1] != 3 or angles.shape[0] == 0:
        raise ValueError("angles must have shape (N,) or (N, 3), got %s" % (angles.shape,))
    n = angles.shape[0]

    dso = _per_view(geo, "DSO", n, ())
    dsd = _per_view(geo, "DSD", n, ())
    off_origin = _per_view(geo, "offOrigin", n, (3,), default=np.zeros(3))  # (z, y, x)
    off_detector = _per_view(geo, "offDetector", n, (2,), default=np.zeros(2))  # (v, u)
    rot_detector = _per_view(geo, "rotDetector", n, (3,), default=np.zeros(3))  # about (x, y, z)
    cor = _per_view(geo, "COR", n, (), default=0.0)
    d_detector = np.asarray(geo.dDetector, dtype=np.float64).reshape(2)  # (v, u)
    n_detector = np.asarray(geo.nDetector, dtype=np.float64).reshape(2)

    rot = _rot_zyz(angles[:, 0], angles[:, 1], angles[:, 2])
    rot_det = _rot_roll_pitch_yaw(rot_detector[:, 2], rot_detector[:, 1], rot_detector[:, 0])
    full = rot @ rot_det

    # Shift of everything relative to the image: COR acts in the rotated x-y
    # plane, offOrigin moves the image (so everything else moves the other way).
    shift = np.stack(
        [-cor * np.sin(angles[:, 0]), cor * np.cos(angles[:, 0]), np.zeros(n)], axis=1
    ) - off_origin[:, ::-1]

    source_local = np.stack([dso, np.zeros(n), np.zeros(n)], axis=1)
    centre_local = np.stack([dso - dsd, off_detector[:, 1], off_detector[:, 0]], axis=1)
    source = np.einsum("nij,nj->ni", rot, source_local) + shift
    centre = np.einsum("nij,nj->ni", rot, centre_local) + shift
    u = d_detector[1] * full[:, :, 1]
    v = d_detector[0] * full[:, :, 2]

    if detector_reference == "first_pixel":
        centre = centre - u * (n_detector[1] - 1) / 2 - v * (n_detector[0] - 1) / 2
    return np.concatenate([source, centre, u, v], axis=1)


def from_vector_to_TIGRE(geo, vectors, detector_reference="center", rtol=1e-6):
    """TIGRE cone-beam geometry from per-view source and detector vectors.

    Parameters
    ----------
    geo : Geometry
        Supplies the image (``nVoxel``, ``sVoxel``, ``dVoxel``) and
        ``nDetector``. It is not modified.
    vectors : numpy.ndarray
        Array of shape ``(N, 12)``. See the module documentation for the
        layout.
    detector_reference : {"center", "first_pixel"}, optional
        Whether columns 3-5 hold the centre of the detector (default, as in
        ASTRA) or the centre of pixel ``(row 0, column 0)``.
    rtol : float, optional
        Relative tolerance for the checks that ``u`` and ``v`` are
        perpendicular and that their lengths are the same in every view.
        Default ``1e-6``.

    Returns
    -------
    Geometry
        A copy of ``geo`` with ``mode="cone"`` and per-view ``DSO``, ``DSD``,
        ``offDetector``, ``rotDetector``, ``offOrigin`` (zero), ``COR`` (zero)
        and ``angles`` of shape ``(N, 3)``, plus ``dDetector`` and
        ``sDetector`` from the lengths of ``u`` and ``v``. Pass
        ``geo.angles`` as the angles to the operators and algorithms.

    Raises
    ------
    ValueError
        If the vectors have the wrong shape or are not finite, if ``u`` and
        ``v`` are not perpendicular, if their lengths change between views,
        if the source is at the centre of the image or in the plane of the
        detector.

    See Also
    --------
    from_TIGRE_to_vector
    """
    _check_reference(detector_reference)
    vectors = np.asarray(vectors, dtype=np.float64)
    if vectors.ndim == 1 and vectors.shape[0] == 12:
        vectors = vectors[None, :]
    if vectors.ndim != 2 or vectors.shape[1] != 12 or vectors.shape[0] == 0:
        raise ValueError("vectors must have shape (N, 12), got %s" % (vectors.shape,))
    if not np.isfinite(vectors).all():
        raise ValueError("vectors must contain only finite values")
    n = vectors.shape[0]
    source, centre, u, v = vectors[:, 0:3], vectors[:, 3:6], vectors[:, 6:9], vectors[:, 9:12]

    n_detector = np.asarray(geo.nDetector, dtype=np.float64).reshape(2)
    if detector_reference == "first_pixel":
        centre = centre + u * (n_detector[1] - 1) / 2 + v * (n_detector[0] - 1) / 2

    # Pixel size: TIGRE has one (dV, dU) for the whole scan, and square axes.
    du = np.linalg.norm(u, axis=1)
    dv = np.linalg.norm(v, axis=1)
    if np.any(du == 0) or np.any(dv == 0):
        raise ValueError("u and v must be non-zero in every view")
    for name, d in (("u", du), ("v", dv)):
        if np.max(np.abs(d - d[0])) > rtol * d[0]:
            raise ValueError(
                "The length of %s (the pixel size) changes between views, from %g to %g. "
                "TIGRE uses one detector pixel size for the whole scan." % (name, d.min(), d.max())
            )
    u_hat = u / du[:, None]
    v_hat = v / dv[:, None]
    cos_uv = np.abs(np.sum(u_hat * v_hat, axis=1))
    if np.any(cos_uv > rtol):
        i = int(np.argmax(cos_uv))
        raise ValueError(
            "u and v are not perpendicular in view %d (angle %.6f degrees). "
            "TIGRE detectors have perpendicular axes."
            % (i, np.degrees(np.arccos(np.clip(np.sum(u_hat[i] * v_hat[i]), -1, 1))))
        )
    # Remove the round-off part of v along u, so the detector frame is exactly a rotation.
    v_hat = v_hat - np.sum(u_hat * v_hat, axis=1)[:, None] * u_hat
    v_hat /= np.linalg.norm(v_hat, axis=1)[:, None]

    dso = np.linalg.norm(source, axis=1)
    if np.any(dso <= rtol * np.max(np.abs(vectors[:, 0:6]))):
        raise ValueError("The source must not be at the centre of the image")
    normal = np.cross(u_hat, v_hat)
    if np.any(np.abs(np.sum((source - centre) * normal, axis=1)) <= rtol * dso):
        raise ValueError("The source must not be in the plane of the detector")

    # Rotation of the TIGRE "angles" frame: its x axis points at the source. The
    # rotation about that axis is free; pick the one whose z axis is closest to
    # world z, so a source in the x-y plane gives angles (alpha, 0, 0).
    x_axis = source / dso[:, None]
    z_axis = np.zeros_like(x_axis)
    for i in range(n):
        for candidate in (np.array([0.0, 0.0, 1.0]), v_hat[i], u_hat[i]):
            perp = candidate - np.dot(candidate, x_axis[i]) * x_axis[i]
            norm = np.linalg.norm(perp)
            if norm > 1e-3:
                z_axis[i] = perp / norm
                break
    y_axis = np.cross(z_axis, x_axis)
    rot = np.stack([x_axis, y_axis, z_axis], axis=2)  # columns are the axes

    angles = _euler_zyz(rot)
    rot = _rot_zyz(angles[:, 0], angles[:, 1], angles[:, 2])  # exact rotation from the angles

    centre_local = np.einsum("nji,nj->ni", rot, centre)  # rot^T @ centre
    detector_frame = np.stack([normal, u_hat, v_hat], axis=2)
    rot_det = np.transpose(rot, (0, 2, 1)) @ detector_frame
    roll, pitch, yaw = _roll_pitch_yaw(rot_det)

    ngeo = copy.deepcopy(geo)
    ngeo.mode = "cone"
    ngeo.DSO = dso
    ngeo.DSD = dso - centre_local[:, 0]
    ngeo.offDetector = np.stack([centre_local[:, 2], centre_local[:, 1]], axis=1)  # (v, u)
    ngeo.rotDetector = np.stack([yaw, pitch, roll], axis=1)  # about (x, y, z)
    ngeo.offOrigin = np.zeros((n, 3))
    ngeo.COR = np.zeros(n)
    ngeo.dDetector = np.array([dv[0], du[0]])
    ngeo.sDetector = ngeo.nDetector * ngeo.dDetector
    ngeo.angles = angles
    return ngeo


def _check_reference(detector_reference):
    if detector_reference not in _REFERENCES:
        raise ValueError(
            "detector_reference must be one of %s, got %r" % (_REFERENCES, detector_reference)
        )


def _per_view(geo, name, n, shape, default=None):
    """``geo.<name>`` as float64 of shape ``(n,) + shape``, given once or per view."""
    value = getattr(geo, name, None)
    if value is None:
        if default is None:
            raise ValueError("geo.%s is not set" % name)
        value = default
    value = np.asarray(value, dtype=np.float64)
    if value.shape in (shape, (1,) + shape) or (shape == () and value.size == 1):
        return np.broadcast_to(value.reshape(shape), (n,) + shape).copy()
    if value.shape == (n,) + shape:
        return value.copy()
    if shape == () and value.size == n:  # e.g. (n, 1)
        return value.reshape(n).copy()
    raise ValueError(
        "geo.%s must have shape %s or %s, got %s" % (name, shape, (n,) + shape, value.shape)
    )


def _rz(a):
    c, s, z, o = np.cos(a), np.sin(a), np.zeros_like(a), np.ones_like(a)
    return np.stack([np.stack([c, -s, z], -1), np.stack([s, c, z], -1), np.stack([z, z, o], -1)], -2)


def _ry(a):
    c, s, z, o = np.cos(a), np.sin(a), np.zeros_like(a), np.ones_like(a)
    return np.stack([np.stack([c, z, s], -1), np.stack([z, o, z], -1), np.stack([-s, z, c], -1)], -2)


def _rx(a):
    c, s, z, o = np.cos(a), np.sin(a), np.zeros_like(a), np.ones_like(a)
    return np.stack([np.stack([o, z, z], -1), np.stack([z, c, -s], -1), np.stack([z, s, c], -1)], -2)


def _rot_zyz(alpha, theta, psi):
    """The kernels' ``eulerZYZ``: Rz(alpha) Ry(theta) Rz(psi)."""
    return _rz(alpha) @ _ry(theta) @ _rz(psi)


def _rot_roll_pitch_yaw(roll, pitch, yaw):
    """The kernels' ``rollPitchYaw``: Rz(roll) Ry(pitch) Rx(yaw)."""
    return _rz(roll) @ _ry(pitch) @ _rx(yaw)


def _euler_zyz(rot, eps=1e-9):
    """Angles (alpha, theta, psi) with Rz(alpha) Ry(theta) Rz(psi) == rot."""
    theta = np.arccos(np.clip(rot[:, 2, 2], -1.0, 1.0))
    sin_theta = np.sin(theta)
    general = sin_theta > eps
    alpha = np.where(
        general,
        np.arctan2(rot[:, 1, 2], rot[:, 0, 2]),
        # theta = 0: Rz(alpha + psi); theta = pi: Rz(alpha) diag(-1, 1, -1) Rz(psi)
        np.where(rot[:, 2, 2] > 0, np.arctan2(rot[:, 1, 0], rot[:, 0, 0]),
                 np.arctan2(-rot[:, 1, 0], rot[:, 1, 1])),
    )
    psi = np.where(general, np.arctan2(rot[:, 2, 1], -rot[:, 2, 0]), 0.0)
    return np.stack([alpha, theta, psi], axis=1)


def _roll_pitch_yaw(rot, eps=1e-9):
    """Angles with Rz(roll) Ry(pitch) Rx(yaw) == rot."""
    pitch = -np.arcsin(np.clip(rot[:, 2, 0], -1.0, 1.0))
    general = np.cos(pitch) > eps
    yaw = np.where(general, np.arctan2(rot[:, 2, 1], rot[:, 2, 2]), 0.0)
    roll = np.where(
        general, np.arctan2(rot[:, 1, 0], rot[:, 0, 0]), np.arctan2(-rot[:, 0, 1], rot[:, 1, 1])
    )
    return roll, pitch, yaw
