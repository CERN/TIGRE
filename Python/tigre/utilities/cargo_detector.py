"""Source-facing cargo DR modules with ray-aligned neighbouring edges."""
import numpy as np


def cargo_module_layout(clearance_mm=20.0, dso_mm=4230.0, horizontal_lift_mm=300.0):
    """Return world centres, pitch angles and arm counts for 36x4, 5 mm modules.

    Lower/upper long-axis edges share source-ray planes over their common
    transverse coverage. Full modules slightly extend the nominal arm envelope.
    DSO defaults to 4230 mm. Starting from the translated reference support,
    fit half a module below the container silhouette, extend the upper arm
    to cover it, then add one full air module. Refit modules face the source.
    The horizontal centre support is lifted 300 mm from its original plane;
    its rearward position and module tilts are refitted to the same joint ray.
    """
    if not np.isfinite(clearance_mm) or clearance_mm < 0:
        raise ValueError("Cargo arm clearance must be finite and nonnegative.")
    if not np.isfinite(dso_mm) or dso_mm <= 0:
        raise ValueError("Cargo DSO must be positive and finite.")
    if not np.isfinite(horizontal_lift_mm) or horizontal_lift_mm < 0:
        raise ValueError("Horizontal arm lift must be finite and nonnegative.")
    source = np.array([4700*np.cos(np.deg2rad(13)), 0., -4700*np.sin(np.deg2rad(13))])
    x = -1800*np.cos(np.deg2rad(39))
    half = 36*5/2
    shift = source*(dso_mm/4700.0-1.0)
    corners = np.array([[cx, 0., cz] for cx in (-1220., 1220.) for cz in (-1195., 1395.)])-shift
    def fan_angle(point):
        ray = point-source
        return np.arctan2(ray[2], -ray[0])
    low, high = min(map(fan_angle, corners)), max(map(fan_angle, corners))

    def pose(center):
        normal = source-center
        normal /= np.linalg.norm(normal)
        row = np.cross(normal, [0., 1., 0.])
        return row, np.rad2deg(np.arctan2(-normal[2], normal[0]))

    def next_module(edge, axis, fixed):
        ray = edge-source
        def evaluate(value):
            center = np.array([fixed, 0., value]) if axis == 2 else np.array([value, 0., fixed])
            row, angle = pose(center)
            lower = center-half*row-source
            residual = lower[0]*ray[2]-lower[2]*ray[0]
            return residual, center, row, angle
        lo, hi = ((-10000., 10000.) if axis == 2 else (-10000., source[0]-1.))
        sign = evaluate(lo)[0]
        if sign*evaluate(hi)[0] >= 0:
            raise ValueError("Cannot align cargo module edge with source ray.")
        for _ in range(70):
            middle = (lo+hi)/2
            if evaluate(middle)[0]*sign > 0:
                lo = middle
            else:
                hi = middle
        _, center, row, angle = evaluate((lo+hi)/2)
        return center, angle, center+half*row

    centers, angles = [], []
    # Aim the centre of the lowest module at the limiting lower corner:
    # its lower half (18 rows) is entirely air in the transverse projection.
    center = np.array([x, 0., source[2]+(source[0]-x)*np.tan(low)])
    row, angle = pose(center)
    centers.append(center); angles.append(angle)
    edge = center+half*row
    vertical_count = round((1024-410)*5/(36*5))+1
    for _ in range(vertical_count-1):
        center, angle, edge = next_module(edge, 2, x)
        centers.append(center); angles.append(angle)
    # If lowering makes the vertical support too short, extend only its top.
    while edge[2]+shift[2] < 1395.0:
        if vertical_count >= 128:
            raise ValueError("Cargo vertical arm exceeds module limit.")
        center, angle, edge = next_module(edge, 2, x)
        centers.append(center); angles.append(angle); vertical_count += 1
    joint = edge.copy()
    near_height = joint[2]+clearance_mm
    edge = source+(joint-source)*((near_height-source[2])/(joint[2]-source[2]))
    # The first lower edge is exactly at the raised, rearward joint point.
    # Its centre lies half a module farther along the source-facing tangent.
    ray = edge-source
    length = np.linalg.norm(ray)
    theta = np.arctan2(ray[2], ray[0])-np.arcsin(half/length)
    radius = np.sqrt(length*length-half*half)
    center = source+radius*np.array([np.cos(theta), 0., np.sin(theta)])
    row, angle = pose(center)
    horizontal_height = center[2]+horizontal_lift_mm
    if horizontal_lift_mm:
        center, angle, edge = next_module(joint, 0, horizontal_height)
    else:
        edge = center+half*row
    centers.append(center); angles.append(angle)
    target_x = x+410*5
    horizontal_count = 1
    while edge[0] < target_x:
        if horizontal_count >= 128:
            raise ValueError("Cargo arm exceeds module limit.")
        center, angle, edge = next_module(edge, 0, horizontal_height)
        centers.append(center); angles.append(angle); horizontal_count += 1
    # Cover the complete transverse box at the requested DSO, without moving
    # existing modules. Conveyor motion supplies the longitudinal scan.
    while fan_angle(edge) < high:
        center, angle, edge = next_module(edge, 0, horizontal_height)
        centers.append(center); angles.append(angle); horizontal_count += 1
        if horizontal_count >= 127:
            raise ValueError("Cargo horizontal arm exceeds module limit.")
    # One complete extra module lies outside the transverse silhouette.
    center, angle, edge = next_module(edge, 0, horizontal_height)
    centers.append(center); angles.append(angle); horizontal_count += 1

    # Express the fitted support in the world frame of the requested source.
    return np.asarray(centers)+shift, np.asarray(angles), vertical_count, horizontal_count
