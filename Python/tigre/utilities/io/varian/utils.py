import numpy as np
from scipy.interpolate import griddata


def get_xmlns(xml_root):
    """Get the namespace of xml file. Utility function for xml readers."""
    return {"": xml_root.tag.split("}")[0].strip("{")}


def _interp_weight(x, xp, N=360):
    """Linear interpolation of xp at x, mod(N). xp must be monotonically increasing"""

    i_min = np.argmin(abs(x - xp))
    if x - xp[i_min] >= 0:
        i_lower = i_min
    else:
        i_lower = (i_min - 1) % len(xp)

    i_upper = (i_lower + 1) % len(xp)

    xp_diff = (xp[i_upper] - xp[i_lower]) % N
    weights = (np.array([(xp[i_upper] - x), (x - xp[i_lower])]) % N) / float(xp_diff)
    return i_lower, i_upper, weights


def interpolate_blank_scan(angle, blank_projs, blank_angles, blank_airnorms):
    """Interpolate blank projections and airnorm values at a given angle. Used in Varian v2.7"""
    i_lower, i_upper, w = _interp_weight(angle, blank_angles)
    blank_interp = w[0] * blank_projs[i_lower] + w[1] * blank_projs[i_upper]
    airnorm_interp = w[0] * blank_airnorms[i_lower] + w[1] * blank_airnorms[i_upper]
    return blank_interp, airnorm_interp


def get_edge_mask(a):
    "Returns an edge mask. Edge values of an N-dim array are set to true."
    a = np.squeeze(a)
    indices = np.indices(a.shape)
    mask = np.zeros(a.shape, dtype=bool)
    for i in range(len(a.shape)):
        edge_indices = [0, a.shape[i] - 1]
        for j in edge_indices:
            mask[indices[i] == j] = 1
    return mask


def mask_outside(a, min_val=0.0, max_val=1.0, clip_edges=True, mask_invalid=True):
    a = np.array(a, copy=True)
    if clip_edges:
        edge_mask = get_edge_mask(a)
        a[edge_mask] = np.clip(a[edge_mask], min_val, max_val)
    if mask_invalid:
        mask = np.ma.masked_outside(a, min_val, max_val).mask | np.ma.masked_invalid(a).mask
        return np.ma.masked_array(a, mask=mask)
    else:
        return np.ma.masked_outside(a, min_val, max_val)


def interp_masked_array(masked_array, fill_value=0.0):
    """Linearly interpolates interior values of an array given a mask. Boundary values assigned
    fill_value (default: 0)."""
    a = np.ma.getdata(masked_array)
    mask = masked_array.mask
    a_interp = a.copy()
    u, v = np.indices(a.shape)
    points = (u[~mask], v[~mask])
    xi = (u[mask], v[mask])
    a_interp[mask] = griddata(points, a[~mask], xi, method="linear", fill_value=fill_value)
    return a_interp
