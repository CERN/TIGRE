import numpy as np
import copy
from types import SimpleNamespace
from scipy.spatial.transform import Rotation
from .geometry import Geometry
from .vector_geometry import from_TIGRE_to_vector, from_vector_to_TIGRE

# Tomosymthesis
def staticDetectorGeo(geo,angles,rot=0) -> Geometry:
    """Geometry for a fixed detector and a source that moves on an arc.

    As in digital breast tomosynthesis: the detector stays still and the
    source moves around the origin. Per-angle ``DSD``, ``offDetector`` and
    ``rotDetector`` are set to describe it.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry for the central position (it is not modified).
    angles : numpy.ndarray
        Angle of the source, in radians, measured from the line through the
        detector centre perpendicular to the detector. 0 means the source
        faces the detector.
    rot : float, optional
        Rotation of the whole setup around the origin, in degrees. Default 0.

    Returns
    -------
    Geometry
        A new geometry. Pass ``geo.angles`` (shape ``(N, 3)``) as the angles to
        the operators and algorithms.

    See Also
    --------
    staticDetLinearSourceGeo, ArbitrarySourceDetMoveGeo
    """
    ngeo = copy.deepcopy(geo)
    R = ngeo.DSD-ngeo.DSO
    rot = rot/180*np.pi
    ngeo.DSD = ngeo.DSO + R*np.cos(angles)
    ngeo.offDetector = np.vstack([0*angles, R*np.sin(angles)]).T
    ngeo.rotDetector = np.vstack([0*angles, 0*angles, -angles]).T
    angles = angles + rot  # do not modify the caller's array
    ngeo.angles = np.vstack([angles, 0*angles, 0*angles]).T
    
    return ngeo

# Linear scan of source
def staticDetLinearSourceGeo(geo,s_pos,s_rot=0,rot=0) -> Geometry:
    """Geometry for a fixed detector and a source that moves along a line.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry for the central position (it is not modified).
        Its detector, including ``offDetector`` and ``rotDetector``, stays
        where it is at angle 0 for every source position.
    s_pos : numpy.ndarray
        Positions of the source along its line, in mm. At 0 the source, origin
        and detector centre are aligned on the x axis.
    s_rot : float or numpy.ndarray, optional
        Angle between the source's line and the detector, in degrees, seen
        from above (x-y plane), anticlockwise. Default 0.
    rot : float or numpy.ndarray, optional
        Rotation of source and detector around the origin, in degrees, seen
        from above, anticlockwise. Default 0.

    Returns
    -------
    Geometry
        A new geometry. Pass ``geo.angles`` (shape ``(N, 3)``) as the angles to
        the operators and algorithms.

    See Also
    --------
    staticDetectorGeo, ArbitrarySourceDetMoveGeo, from_vector_to_TIGRE
    """
    s_pos = np.atleast_1d(np.asarray(s_pos, dtype=np.float64))
    if s_pos.ndim != 1:
        raise ValueError("s_pos should be a 1D array of positions along the source's line")
    n = s_pos.shape[0]
    s_rot = np.radians(_broadcast_per_view(s_rot, n, "s_rot"))
    rot = np.radians(_broadcast_per_view(rot, n, "rot"))

    central = _central_vectors(geo)
    vectors = np.tile(central, (n, 1))
    vectors[:, 0] += s_pos * np.sin(s_rot)
    vectors[:, 1] += s_pos * np.cos(s_rot)
    return from_vector_to_TIGRE(geo, _rotate_about_origin(vectors, rot, _central_offset(geo)))


def ArbitrarySourceDetMoveGeo(geo,s_pos,d_pos=None,d_rot=None) -> Geometry:
    """Geometry for a source and detector that move freely around a fixed object.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry (it is not modified). It gives the image and the
        detector's pixels; ``offOrigin`` places the image as usual.
    s_pos : numpy.ndarray
        Source positions ``(x, y, z)`` in mm, shape ``(N, 3)``. The source of
        ``geo`` sits at ``(DSO, 0, 0)``.
    d_pos : numpy.ndarray, optional
        Detector centre positions ``(x, y, z)`` in mm, shape ``(N, 3)``.
        Default ``(DSO - DSD, 0, 0)`` for every position.
    d_rot : numpy.ndarray, optional
        Detector rotations ``(roll, pitch, yaw)`` in degrees, shape ``(N, 3)``,
        in the same order as ``geo.rotDetector``, to which they are added
        (roll turns the detector within its own plane). Default: no rotation,
        so the detector faces the origin, with its rows as horizontal as
        possible.

    Returns
    -------
    Geometry
        A new geometry. Pass ``geo.angles`` (shape ``(N, 3)``) as the angles to
        the operators and algorithms.

    See Also
    --------
    staticDetectorGeo, staticDetLinearSourceGeo, from_vector_to_TIGRE
    """
    s_pos = np.asarray(s_pos, dtype=np.float64)
    if s_pos.ndim != 2 or s_pos.shape[1] != 3 or s_pos.shape[0] == 0:
        raise ValueError("s_pos should be an N x 3 array, got shape %s" % (s_pos.shape,))
    n = s_pos.shape[0]
    first = _first_view(geo)
    if d_pos is None:
        d_pos = np.tile([first.DSO - first.DSD, 0.0, 0.0], (n, 1))
    d_pos = np.asarray(d_pos, dtype=np.float64)
    d_rot = np.zeros((n, 3)) if d_rot is None else np.asarray(d_rot, dtype=np.float64)
    if d_pos.shape != s_pos.shape or d_rot.shape != s_pos.shape:
        raise ValueError(
            "s_pos, d_pos and d_rot should all be N x 3, got %s, %s and %s"
            % (s_pos.shape, d_pos.shape, d_rot.shape)
        )
    if np.any(np.linalg.norm(d_pos, axis=1) == 0):
        raise ValueError("The detector centre must not be at the origin")

    vectors = np.zeros((n, 12))
    rot_detector = first.rotDetector + np.radians(d_rot)  # (roll, pitch, yaw): about the detector's (x, y, z)
    for i in range(n):
        base = _facing_origin(d_pos[i])
        frame = base @ _rz(rot_detector[i, 2]) @ _ry(rot_detector[i, 1]) @ _rx(rot_detector[i, 0])
        vectors[i, 0:3] = s_pos[i]
        vectors[i, 3:6] = d_pos[i] + base @ np.array([0.0, first.offDetector[1], first.offDetector[0]])
        vectors[i, 6:9] = first.dDetector[1] * frame[:, 1]
        vectors[i, 9:12] = first.dDetector[0] * frame[:, 2]
    vectors[:, 0:6] -= np.tile(first.offOrigin[::-1], 2)  # relative to the image centre
    return from_vector_to_TIGRE(geo, vectors)


def ArbitrarySourceDetectorFixedObject(
        geometry: Geometry,
        focal_spot_position_mm: np.ndarray,
        detector_center_position_mm: np.ndarray,
        detector_line_direction: np.ndarray,
        detector_column_direction: np.ndarray,
        origin_mm: np.ndarray = None,
        use_center_correction: bool = True) -> Geometry:
    """Geometry from per-view source and detector positions in a scanner's own frame.

    Parameters
    ----------
    geometry : Geometry
        Gives the image and the detector's pixels (``nDetector``,
        ``dDetector``). It is updated in place and returned.
    focal_spot_position_mm : numpy.ndarray
        Source positions, shape ``(N, 3)``.
    detector_center_position_mm : numpy.ndarray
        Detector centre positions, shape ``(N, 3)``.
    detector_line_direction : numpy.ndarray
        Direction from pixel ``(0, 0)`` to ``(0, 1)`` (along a row), shape ``(N, 3)``.
    detector_column_direction : numpy.ndarray
        Direction from pixel ``(0, 0)`` to ``(1, 0)`` (down a column), shape ``(N, 3)``.
    origin_mm : numpy.ndarray, optional
        Subtracted from the positions, so the image is centred there. Default 0.
    use_center_correction : bool, optional
        If True (default), also subtract the point that is closest, in the
        least-squares sense, to all the source-detector lines, so the image
        is centred on the trajectory.

    Returns
    -------
    Geometry
        ``geometry``, with per-view ``DSO``, ``DSD``, ``offDetector``,
        ``rotDetector``, ``offOrigin`` (zero), ``COR`` (zero) and ``angles``.

    See Also
    --------
    from_vector_to_TIGRE
    """
    source = np.asarray(focal_spot_position_mm, dtype=np.float64)
    centre = np.asarray(detector_center_position_mm, dtype=np.float64)
    line = np.asarray(detector_line_direction, dtype=np.float64)
    column = np.asarray(detector_column_direction, dtype=np.float64)
    if origin_mm is not None:
        source = source - np.asarray(origin_mm, dtype=np.float64).reshape(1, 3)
        centre = centre - np.asarray(origin_mm, dtype=np.float64).reshape(1, 3)
    if use_center_correction:
        trajectory_center_mm = calculate_trajectory_center_mm(source, centre)
        source = source - trajectory_center_mm
        centre = centre - trajectory_center_mm

    d_detector = np.asarray(geometry.dDetector, dtype=np.float64).reshape(2)  # (v, u)
    u = d_detector[1] * line / np.linalg.norm(line, axis=1, keepdims=True)
    v = d_detector[0] * column / np.linalg.norm(column, axis=1, keepdims=True)
    converted = from_vector_to_TIGRE(geometry, np.concatenate([source, centre, u, v], axis=1))
    for name in ("mode", "DSO", "DSD", "offDetector", "rotDetector", "offOrigin", "COR", "angles"):
        setattr(geometry, name, getattr(converted, name))
    return geometry


def _first_view(geo):
    """The detector and distances of ``geo``'s first (or only) view, as float64 arrays."""
    def value(name, shape, default):
        v = getattr(geo, name, None)
        v = np.asarray(default if v is None else v, dtype=np.float64)
        return v.reshape((-1,) + shape)[0] if v.size > int(np.prod(shape)) else v.reshape(shape)

    return SimpleNamespace(
        mode="cone",
        DSO=float(value("DSO", (), None)),
        DSD=float(value("DSD", (), None)),
        offOrigin=value("offOrigin", (3,), np.zeros(3)),
        offDetector=value("offDetector", (2,), np.zeros(2)),
        rotDetector=value("rotDetector", (3,), np.zeros(3)),
        COR=float(value("COR", (), 0.0)),
        dDetector=np.asarray(geo.dDetector, dtype=np.float64).reshape(2),
        nDetector=np.asarray(geo.nDetector).reshape(2),
    )


def _central_vectors(geo):
    """Vectors of ``geo``'s first view at angle 0."""
    return from_TIGRE_to_vector(_first_view(geo), np.zeros((1, 3)))[0]


def _central_offset(geo):
    """``offOrigin`` of the first view as (x, y, z)."""
    return _first_view(geo).offOrigin[::-1]


def _broadcast_per_view(value, n, name):
    value = np.asarray(value, dtype=np.float64)
    if value.ndim == 0:
        return np.full(n, float(value))
    if value.shape != (n,):
        raise ValueError("%s should be a scalar or have one value per source position" % name)
    return value


def _rotate_about_origin(vectors, angles, image_centre):
    """Rotate scans (given relative to the image centre) about the z axis through the origin."""
    out = np.empty_like(vectors)
    for i, a in enumerate(angles):
        r = _rz(a)
        out[i, 0:3] = r @ (vectors[i, 0:3] + image_centre) - image_centre
        out[i, 3:6] = r @ (vectors[i, 3:6] + image_centre) - image_centre
        out[i, 6:9] = r @ vectors[i, 6:9]
        out[i, 9:12] = r @ vectors[i, 9:12]
    return out


def _facing_origin(position):
    """Rotation taking TIGRE's unrotated detector (normal +x, rows along +y) to one at ``position`` facing the origin."""
    normal = -position / np.linalg.norm(position)
    up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(up, normal)) > 1 - 1e-6:
        up = np.array([1.0, 0.0, 0.0])
    up = up - np.dot(up, normal) * normal
    up /= np.linalg.norm(up)
    return np.stack([normal, np.cross(up, normal), up], axis=1)


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _rx(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def euler_from_vecs(a_vec,b_vec,order="xyz"):
    """
    Calculate Euler angles from two vectors

    Parameters
    ----------
    a_vec : np.array of (n,3) or (3,)
        A vector(s).
    b_vec : np.array of (n,3) or (3,)
        B vector(s).
    order : string, optional
        Order of the euler angles in rotations. The default is "xyz".

    Returns
    -------
    euler : np.array of (n,3) or (3,)
        Euler angles of rotations.

    """
    na = 1 if a_vec.ndim<2 else a_vec.shape[0]
    nb = 1 if b_vec.ndim<2 else b_vec.shape[0]
    n = max(na,nb)
    a_vec = np.repeat([a_vec], n, axis=0) if a_vec.ndim<2 else a_vec
    b_vec = np.repeat([b_vec], n, axis=0) if b_vec.ndim<2 else b_vec
    euler = np.zeros_like(a_vec,dtype=np.float64)
    for i in range(n):
        R = rotation_from_vecs(a_vec[i,:],b_vec[i,:])
        euler[i,:] = Rotation.from_matrix(R).as_euler(order)
        for j in range(3):
            if abs(euler[i,j]) < 2e-8: # very small angle
                euler[i,j] = 0 
    return euler


def rotation_from_vecs(v1, v2):
    """
    Compute a matrix R that rotates v1 to align with v2.
    v1 and v2 must be length-3 1d numpy arrays.
    """
    # unit vectors
    u = v1 / np.linalg.norm(v1)
    Ru = v2 / np.linalg.norm(v2)
    # dimension of the space and identity
    I = np.identity(u.size)
    # the cos angle between the vectors (dot product)
    c = np.dot(u, Ru)
    # a small number
    eps = 1.0e-9
    if np.abs(c - 1.0) < eps:
        # same direction
        return I
    elif np.abs(c + 1.0) < eps:
        # opposite direction: half a turn about any axis perpendicular to u (-I would be a reflection)
        e = I[np.argmin(np.abs(u))]
        a = np.cross(u, e)
        a /= np.linalg.norm(a)
        return 2 * np.outer(a, a) - I
    else:
        # the cross product matrix of a vector to rotate around
        v = np.cross(u, Ru)
        K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        # Rodrigues' formula
        return I + K + (K @ K) / (1 + c)    


def perpendicular_point_on_line(point: np.ndarray, line_point: np.ndarray, line_direction: np.ndarray):
    ap = point - line_point
    dot_product = np.dot(ap, line_direction)
    line_squared = np.dot(line_direction, line_direction)
    return line_point + (dot_product / line_squared) * line_direction

def calculate_trajectory_center_mm(focal_spot_position_mm: np.ndarray, detector_center_position_mm: np.ndarray):
    number_of_projection = focal_spot_position_mm.shape[0]
    A = np.zeros((number_of_projection * 6, 5))
    b = np.zeros((number_of_projection * 6, 1))


    for i in range(0, number_of_projection * 6-1, 6):
        ii = i // 6
        source = focal_spot_position_mm[ii]
        detector = detector_center_position_mm[ii]
        direction = source - detector
        direction = direction / np.linalg.norm(direction)

        A[i, :] = np.array([1, 0, 0, direction[0], 0])
        b[i] = source[0]
        i=i+1

        A[i, :] = np.array([0, 1, 0, direction[1], 0])
        b[i] = source[1]
        i=i+1

        A[i, :] = np.array([0, 0, 1, direction[2], 0])
        b[i] = source[2]
        i=i+1

        A[i, :] = np.array([1, 0, 0, 0, -direction[0]])
        b[i] = detector[0]
        i=i+1

        A[i, :] = np.array([0, 1, 0, 0, -direction[1]])
        b[i] = detector[1]
        i=i+1

        A[i, :] = np.array([0, 0, 1, 0, -direction[2]])
        b[i] = detector[2]
        i=i+1            
    
    res = np.linalg.lstsq(A, b, rcond=None)
    trajectory_center_mm = np.array([res[0][0], res[0][1], res[0][2]]).reshape((1, 3))
    return trajectory_center_mm

