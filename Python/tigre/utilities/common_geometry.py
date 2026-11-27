import numpy as np
import copy
from scipy.spatial.transform import Rotation
from .geometry import Geometry

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

def tiltedAxisGeo(geo, cor_top, cor_bottom, height=None) -> Geometry:
    """Geometry for a rotation axis tilted within the detector plane.

    Describes a rotation axis that is not parallel to the detector columns:
    its centre of rotation is ``cor_top`` at the top of the image volume (+z)
    and ``cor_bottom`` at the bottom (-z), for example two values of
    :func:`tigre.computeCOR` at two detector rows. The top is the end seen by
    the last detector rows, ``proj[:, -1, :]``; row 0 is the bottom.

    The tilt is described exactly with fields TIGRE already has: ``COR`` is
    replaced, the tilt angle is added to the in-plane detector rotation
    ``rotDetector[..., 0]``, and a small shift is added to the image offset
    along the rotation axis, ``offOrigin[..., 0]``. The reconstructed image is
    aligned with the rotation axis.

    Iterative algorithms are exact. FDK filters along the detector rows, so it
    is approximate for a tilted axis; for large tilts, rotate the projections
    by the tilt angle instead and set only ``COR``.

    Parameters
    ----------
    geo : Geometry
        The scanner geometry (it is not modified).
    cor_top, cor_bottom : float
        Centre of rotation at the top and at the bottom, in mm at the rotation
        axis, with the sign that ``geo.COR`` uses.
    height : float, optional
        Distance between the two measurements along the rotation axis, in mm.
        Default ``geo.sVoxel[0]``, the height of the image volume.

    Returns
    -------
    Geometry
        A copy of ``geo`` with ``COR``, ``rotDetector`` and ``offOrigin`` set.

    Examples
    --------
    >>> geo = tigre.tiltedAxisGeo(geo, cor_top=1.2, cor_bottom=-0.4)
    """
    geo = copy.deepcopy(geo)
    if height is None:
        height = float(np.asarray(geo.sVoxel)[0])
    if not (np.isscalar(cor_top) and np.isscalar(cor_bottom) and np.isscalar(height)):
        raise ValueError("cor_top, cor_bottom and height should be scalars")
    if height <= 0:
        raise ValueError("height should be positive")
    s = (cor_bottom - cor_top) / height
    if abs(s) >= 1:
        raise ValueError("cor_top and cor_bottom differ by more than height")

    phi = np.arcsin(s)  # tilt of the axis within the detector plane
    cor_mid = (cor_top + cor_bottom) / 2  # centre of rotation half way between the two

    geo.COR = cor_mid * np.cos(phi)
    rot = getattr(geo, "rotDetector", None)
    rot = np.zeros(3) if rot is None else np.array(rot, dtype=np.float64)
    rot[..., 0] += phi
    geo.rotDetector = rot
    off = getattr(geo, "offOrigin", None)
    off = np.zeros(3) if off is None else np.array(off, dtype=np.float64)
    off[..., 0] -= cor_mid * np.sin(phi)  # Python order is (z, y, x)
    geo.offOrigin = off
    return geo


def computeImageGeo(geo) -> Geometry:
    """Image size and voxel size that suit the detector and the scan.

    The voxels get the size of a detector pixel seen at the rotation axis,
    ``dDetector * DSO / DSD`` (``dDetector`` for parallel beam), and there are
    enough of them to cover the field of view: the circle that every
    projection sees in the axial plane (larger with an offset, half-fan,
    detector), and the detector's height at the axis. Useful for real data,
    where the image size is not given.

    Parameters
    ----------
    geo : Geometry
        The scanner geometry, with the detector and distances set (it is not
        modified). Per-projection ``DSO``, ``DSD`` and ``offDetector`` are
        averaged / their largest offset is used.

    Returns
    -------
    Geometry
        A copy of ``geo`` with ``nVoxel``, ``dVoxel`` and ``sVoxel`` set.

    Examples
    --------
    >>> geo = tigre.computeImageGeo(geo)
    """
    geo = copy.deepcopy(geo)
    d_v, d_u = (float(x) for x in np.asarray(geo.dDetector, dtype=np.float64).ravel()[:2])
    n_v, n_u = (float(x) for x in np.asarray(geo.nDetector, dtype=np.float64).ravel()[:2])
    off = np.zeros(2) if getattr(geo, "offDetector", None) is None else np.asarray(geo.offDetector, dtype=np.float64).reshape(-1, 2)
    off_v, off_u = np.abs(off).max(axis=0)
    u_max = n_u * d_u / 2 + off_u  # furthest detector column from the central ray, mm
    v_max = n_v * d_v / 2 + off_v
    if getattr(geo, "mode", "cone") == "parallel":
        mag = 1.0
        radius = u_max
    else:
        DSO = float(np.mean(geo.DSO))
        DSD = float(np.mean(geo.DSD))
        mag = DSD / DSO
        radius = DSO * u_max / np.hypot(DSD, u_max)  # distance of the outermost ray from the axis
    d_xy = d_u / mag
    d_z = d_v / mag
    n_xy = int(np.ceil(2 * radius / d_xy))
    n_z = int(np.ceil(2 * v_max / mag / d_z))
    geo.nVoxel = np.array([n_z, n_xy, n_xy])
    geo.dVoxel = np.array([d_z, d_xy, d_xy])
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    return geo


def changeResolution(geo, factor, proj=None):
    """Coarser version of a geometry, and optionally of its projections.

    The image keeps its physical size, ``sVoxel``, with ``factor`` times fewer
    voxels along each axis. If ``proj`` is given, the projections are binned by
    the same factor (the mean of each ``factor`` x ``factor`` block of pixels)
    and the detector changes to match; pixels that do not fill a whole block,
    at the end of each detector row or column, are dropped, and ``offDetector``
    moves so that the binned pixels stay where they were. Useful to try
    parameters on a small problem first.

    Parameters
    ----------
    geo : Geometry
        The scanner geometry (it is not modified).
    factor : int
        Reduction factor, 1 or more.
    proj : numpy.ndarray, optional
        Projections, shape ``(N, nDetector[0], nDetector[1])``.

    Returns
    -------
    Geometry, or (Geometry, numpy.ndarray)
        The new geometry; with ``proj``, also the binned projections, float32.

    Examples
    --------
    >>> geo_small, proj_small = tigre.changeResolution(geo, 2, proj)
    """
    factor = int(factor)
    if factor < 1:
        raise ValueError("factor should be an integer, 1 or more")
    geo = copy.deepcopy(geo)
    n_voxel = np.asarray(geo.nVoxel)
    geo.nVoxel = np.maximum(np.round(n_voxel / factor), 1).astype(int)
    geo.dVoxel = np.asarray(geo.sVoxel, dtype=np.float64) / geo.nVoxel
    if proj is None:
        return geo

    proj = np.asarray(proj)
    n_det = np.asarray(geo.nDetector).astype(int)
    if proj.ndim != 3 or tuple(proj.shape[1:]) != tuple(n_det):
        raise ValueError("proj should have shape (N, nDetector[0], nDetector[1])")
    n_new = n_det // factor
    if (n_new < 1).any():
        raise ValueError("factor is larger than the detector")
    kept = n_new * factor
    binned = proj[:, : kept[0], : kept[1]].reshape(proj.shape[0], n_new[0], factor, n_new[1], factor)
    binned = binned.mean(axis=(2, 4), dtype=np.float64).astype(np.float32)
    d_det = np.asarray(geo.dDetector, dtype=np.float64)
    # The dropped pixels are at the high-index end; the kept ones are centred lower by half of them.
    # u grows with the column and v with the row (row 0 is the -z end of the detector).
    shift = (n_det - kept) * d_det / 2
    if getattr(geo, "offDetector", None) is not None:
        geo.offDetector = np.asarray(geo.offDetector, dtype=np.float64) - shift
    else:
        geo.offDetector = -shift
    geo.nDetector = n_new
    geo.dDetector = d_det * factor
    geo.sDetector = geo.nDetector * geo.dDetector
    return geo, binned


def spotDriftGeo(geo, drift, angles) -> Geometry:
    """Geometry for a source focal spot that drifts during the scan.

    The drift is the movement of the source from its nominal position, in mm,
    measured in the frame that rotates with the scanner: along the detector
    ``u`` and ``v`` directions, and optionally along the central ray. It is
    described exactly with fields TIGRE already has: the source and the
    detector are moved back together, so the source returns to its nominal
    place, the detector offset ``offDetector`` takes the opposite shift, and
    the image offset ``offOrigin`` takes it too, rotated with each angle. A
    drift along the central ray changes ``DSO`` and ``DSD``. All of these
    become per-projection values.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry (it is not modified).
    drift : numpy.ndarray
        Shape ``(N, 2)`` or ``(N, 3)``, one row per projection, or a single
        row for a constant offset. Columns, in mm: drift along ``v`` (same
        direction as ``offDetector[..., 0]``, +z for a standard scan), along
        ``u`` (same direction as ``offDetector[..., 1]``) and, optionally,
        along the central ray towards the detector.
    angles : numpy.ndarray
        The projection angles, shape ``(N,)`` or ``(N, 3)``, as given to the
        operators and algorithms.

    Returns
    -------
    Geometry
        A copy of ``geo`` with per-projection ``offOrigin``, ``offDetector``,
        ``DSO`` and ``DSD``.

    Examples
    --------
    >>> geo = tigre.spotDriftGeo(geo, drift, angles)  # drift: (N, 2), mm
    """
    if getattr(geo, "mode", "cone") == "parallel":
        raise ValueError("spotDriftGeo needs a cone-beam geometry; a parallel beam has no source position")
    geo = copy.deepcopy(geo)
    angles = np.asarray(angles, dtype=np.float64)
    if angles.ndim == 1:
        angles = np.stack([angles, 0 * angles, 0 * angles], axis=1)
    if angles.ndim != 2 or angles.shape[1] != 3:
        raise ValueError("angles should have shape (N,) or (N, 3)")
    n = angles.shape[0]
    drift = np.atleast_2d(np.asarray(drift, dtype=np.float64))
    if drift.ndim != 2 or drift.shape[1] not in (2, 3) or drift.shape[0] not in (1, n):
        raise ValueError("drift should have shape (N, 2) or (N, 3), or a single row, with N = number of angles")
    drift = np.broadcast_to(drift, (n, drift.shape[1]))
    d_v, d_u = drift[:, 0], drift[:, 1]
    d_w = drift[:, 2] if drift.shape[1] == 3 else np.zeros(n)

    def per_angle(value, width):
        value = np.zeros(width) if value is None else np.asarray(value, dtype=np.float64)
        return np.array(np.broadcast_to(value, (n, width) if width else (n,)))

    # Moving source and detector by -drift brings the source back; the image moves with them.
    # The rotating frame's u and v axes in world coordinates (columns 2 and 3 of the Euler ZYZ
    # rotation the projectors use).
    a, t, p = angles[:, 0], angles[:, 1], angles[:, 2]
    ax_u = np.stack([-np.cos(a) * np.cos(t) * np.sin(p) - np.sin(a) * np.cos(p),
                     -np.sin(a) * np.cos(t) * np.sin(p) + np.cos(a) * np.cos(p),
                     np.sin(t) * np.sin(p)], axis=1)  # (x, y, z)
    ax_v = np.stack([np.cos(a) * np.sin(t), np.sin(a) * np.sin(t), np.cos(t)], axis=1)
    shift_xyz = d_u[:, None] * ax_u + d_v[:, None] * ax_v
    off_origin = per_angle(getattr(geo, "offOrigin", None), 3)
    geo.offOrigin = off_origin - shift_xyz[:, ::-1]  # Python order is (z, y, x)
    off_det = per_angle(getattr(geo, "offDetector", None), 2)
    geo.offDetector = off_det - np.stack([d_v, d_u], axis=1)
    geo.DSO = per_angle(geo.DSO, 0) - d_w
    geo.DSD = per_angle(geo.DSD, 0) - d_w
    return geo


# Linear scan of source
def staticDetLinearSourceGeo(geo,s_pos,s_rot=0,rot=0) -> Geometry:
    """Geometry for a fixed detector and a source that moves along a line.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry for the central position (it is not modified).
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
    staticDetectorGeo, ArbitrarySourceDetMoveGeo
    """
    ngeo = copy.deepcopy(geo)
    if np.isscalar(s_rot):
        s_rot = s_rot * np.ones_like(s_pos)
    s_rot = s_rot / 180 *np.pi
    if np.isscalar(rot):
        rot = rot * np.ones_like(s_pos)
    if len(s_pos) != len(rot) or len(s_pos) != len(s_rot):
        raise("Inputs length do not match")
    
    rot = rot/180*np.pi 
               
    ang  = np.arctan2(s_pos*np.cos(s_rot), ngeo.DSO + s_pos*np.sin(s_rot)) 
    R = ngeo.DSD - ngeo.DSO
    if hasattr(ngeo,"offDetector"):
        ngeo.offDetector = ngeo.offDetector.astype(np.float32) + np.array([ang*0, R*np.sin(ang)]).T 
    else:
        ngeo.offDetector = np.array([ang*0, R*np.sin(ang)]).T 
    if hasattr(ngeo,"rotDetector"):    
        ngeo.rotDetector += np.array([0*ang, 0*ang, -ang]).T   
    else:
        ngeo.rotDetector = np.array([0*ang, 0*ang, -ang]).T 
    ngeo.DSO = np.sqrt((s_pos*np.cos(s_rot))**2 + (ngeo.DSO + s_pos*np.sin(s_rot))**2)
    ngeo.DSD = np.sqrt((s_pos*np.cos(s_rot))**2 + (ngeo.DSD + s_pos*np.sin(s_rot))**2)
    ang += rot    
    ngeo.angles = np.vstack([ang, 0*ang, 0*ang]).T
    return ngeo


def ArbitrarySourceDetMoveGeo(geo,s_pos,d_pos=None,d_rot=None) -> Geometry:
    """Geometry for a source and detector that move freely around a fixed object.

    Parameters
    ----------
    geo : Geometry
        A cone-beam geometry (it is not modified).
    s_pos : numpy.ndarray
        Source positions ``(x, y, z)`` in mm, shape ``(N, 3)``. The source of
        ``geo`` sits at ``(DSO, 0, 0)``.
    d_pos : numpy.ndarray, optional
        Detector centre positions ``(x, y, z)`` in mm, shape ``(N, 3)``.
        Default ``(DSO - DSD, 0, 0)`` for every position.
    d_rot : numpy.ndarray, optional
        Detector rotations ``(roll, pitch, yaw)`` in degrees, shape ``(N, 3)``.
        Default: no rotation, the detector facing the origin.

    Returns
    -------
    Geometry
        A new geometry. Pass ``geo.angles`` (shape ``(N, 3)``) as the angles to
        the operators and algorithms.

    See Also
    --------
    staticDetectorGeo, staticDetLinearSourceGeo
    """
    s_pos = np.array(s_pos) if isinstance(s_pos,list) else s_pos
    ngeo = copy.deepcopy(geo)
    if s_pos.ndim != 2:
        raise("Input s_pos should be an n x 3 array")
    else:
        n = s_pos.shape[0]
    if s_pos.shape[1] != 3:
        raise("Input s_pos should be an n x 3 array")
    if d_pos is None:
        d_pos = np.repeat([[ngeo.DSO-ngeo.DSD, 0, 0]], n, axis=0)
    elif isinstance(d_pos,list):
        d_pos = np.array(d_pos)
    if d_rot is None:
        d_rot = np.zeros((n,3))
    elif isinstance(d_rot,list):
        d_rot = np.array(d_rot)
    if s_pos.shape != d_pos.shape or s_pos.shape != d_rot.shape:
        raise("Inputs dimensions do not match")        
    
    d_rot = d_rot / 180*np.pi
    
    # source and detector vector lengths and directions
    rs = np.linalg.norm(s_pos, axis=1)
    ns = s_pos / np.tile(rs, (3,1)).T
    rd = np.linalg.norm(d_pos, axis=1)
    nd = d_pos / np.tile(rd, (3,1)).T
    
    # source euler angles (intrinsic) away from x-axis in "ZYZ" order
    s0 = np.zeros_like(s_pos)
    s0[:,0] = s_pos[:,0]
    s_ang = euler_from_vecs(s0, s_pos, order="zyz")   
    
    # check angles between OS and OD
    sd_cross = np.cross(ns, nd)
    # note: A x B = ||A||.||B||sin(theta)
    sd_ang = np.arcsin(sd_cross)
    if (abs(sd_ang) >= np.pi/2).any():  
        raise RuntimeError("Source and detector are not always on different sides (some angle<=90 degree)")        
    
    # detector angles from ray
    d_cross = np.cross(-ns, nd) 
    d_ang = np.arcsin(d_cross)
    pd = d_cross * np.tile(rd, (3,1)).T  
    
    if hasattr(ngeo,"offDetector"):
        ngeo.offDetector = ngeo.offDetector.astype(np.float64) + np.vstack([pd[:,1], -pd[:,2]]).T
    else:
        ngeo.offDetector = np.vstack([pd[:,1], -pd[:,2]]).T
        
    if hasattr(ngeo,"rotDetector"):
        ngeo.rotDetector += d_ang + d_rot
    else:
        ngeo.rotDetector = d_ang + d_rot
    ngeo.DSO = np.linalg.norm(s_pos, axis=1) 
    ngeo.DSD = np.linalg.norm(s_pos-d_pos, axis=1) 
    ngeo.angles = s_ang 
    return ngeo


def ArbitrarySourceDetectorFixedObject(
        geometry: Geometry,
        focal_spot_position_mm: np.ndarray, 
        detector_center_position_mm: np.ndarray, 
        detector_line_direction: np.ndarray, 
        detector_column_direction: np.ndarray,
        origin_mm: np.ndarray = None,
        use_center_correction: bool = True) -> Geometry:
    """
    geo: Geometry object
    focal_spot_position_mm: position of the source, 
    detector_center_position_mm: position of the detector center, 
    detector_line_direction: detector line vector from pixel (0, 0) -> (0, 1), 
    detector_column_direction: detector column vector from pixel (0, 0) -> (1, 0),
    origin_mm: origin of the ct trajectory. The source and detector positions are translated with this value. Defaults to: None.
    use_center_correction: Calculate an arbitrary origin of the trajectory. Defaults to: True.
    """

    # Assumption: CT trajectory has one rotation center.
    number_of_projection = focal_spot_position_mm.shape[0]
    
    if origin_mm is None:
        origin_mm = np.zeros((3,))

    focal_spot_position_mm = focal_spot_position_mm - origin_mm
    detector_center_position_mm = detector_center_position_mm - origin_mm
    
    if use_center_correction:
        trajectory_center_mm = calculate_trajectory_center_mm(
            focal_spot_position_mm, detector_center_position_mm)
        focal_spot_position_mm = focal_spot_position_mm - trajectory_center_mm
        detector_center_position_mm = detector_center_position_mm - trajectory_center_mm
    else:
        trajectory_center_mm = np.zeros((3,))

    if not use_center_correction:
        geometry.offOrigin = trajectory_center_mm.reshape((3, ))

    # source and detector are orthogonal. The angle is rotates the source from the x axis.

    # 1. find nearest point from source detector line to trajectory center
    source_detector_vector = detector_center_position_mm - focal_spot_position_mm
    fdd_mm = np.linalg.norm(source_detector_vector, axis=1).reshape((-1, 1))
    fod_mm = np.zeros_like(fdd_mm)
    source_detector_direction = source_detector_vector / fdd_mm
    
    nearest_point_mm = np.zeros_like(focal_spot_position_mm)
    detector_offsets_mm = np.zeros((number_of_projection, 2))
    
    euler_zyz = np.zeros_like(focal_spot_position_mm)
    euler_xyz = np.zeros_like(focal_spot_position_mm)
    first_rot = np.eye(3)

    for i in range(number_of_projection):
        nearest_point_mm[i] = perpendicular_point_on_line(
            trajectory_center_mm, focal_spot_position_mm[i], -source_detector_direction[i])
        fod_mm[i] = np.linalg.norm(nearest_point_mm[i] - focal_spot_position_mm[i])

        # 2. calculate the angle rotation + offset and check it!
        rotation_matrix = rotation_from_vecs(np.array([1, 0, 0]),-nearest_point_mm[i] +focal_spot_position_mm[i])
        if i == 0:
            first_rot = rotation_matrix.T
            first_rot = np.eye(3)
        angle = rotation_matrix @ first_rot
        if not np.isclose(np.linalg.det(rotation_matrix), 1, 0.01):
            raise ValueError('Rotation matrix must be right handed!')
        
        rotation_inverse = Rotation.from_matrix(rotation_matrix.T)

        offsets = rotation_inverse.apply(nearest_point_mm[i])
        fod_mm[i] += offsets[0]
        detector_offsets_mm[i, 0] = offsets[2] * 2
        detector_offsets_mm[i, 1] = offsets[1] * 2

        # calculate the relative rotation of the detector
        detector_matrix = np.eye(3)
        detector_matrix[:, 1] = detector_line_direction[i]
        detector_matrix[:, 2] = detector_column_direction[i]
        detector_matrix[:, 0] = np.cross(detector_matrix[:, 1], detector_matrix[:, 2])

        if not np.isclose(np.linalg.det(detector_matrix), 1, 0.01):
            raise ValueError('Rotation matrix must be right handed!')
        
        # the rotation is from the angle rotation -> "real" rotation
        relative_rotation_matrix = rotation_matrix.T @ detector_matrix
        # relative_rotation_matrix = detector_matrix @ rotation_matrix.T -> wrong
        relative_rotation = Rotation.from_matrix(relative_rotation_matrix)

        euler_zyz[i] = Rotation.from_matrix(angle).as_euler('zyz', False)
        euler_xyz[i] = relative_rotation.as_euler('xyz', False) 
    
    geometry.DSO = fod_mm
    geometry.DSD = fdd_mm
    geometry.offDetector = detector_offsets_mm
    geometry.rotDetector = euler_xyz
    geometry.angles = euler_zyz
    return geometry

    
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
        # opposite direction
        return -I
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

    print(f'Calculated FOD / ODD')

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
    
    res = np.linalg.lstsq(A, b)
    trajectory_center_mm = np.array([res[0][0], res[0][1], res[0][2]]).reshape((1, 3))
    return trajectory_center_mm

