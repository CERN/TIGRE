import numpy as np


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
