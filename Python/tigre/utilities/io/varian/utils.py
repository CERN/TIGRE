import numpy as np


class XML:
    def __init__(self, filepath, xml_reader):
        self.root = xml_reader(filepath)
        self.ns = self._get_xmlns()

    def _get_xmlns(self):
        """Returns the namespace of the xml file."""
        return {"": self.root.tag.split("}")[0].strip("{")}


def interp_weight(x, xp, N=360):
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


def mm2cm(x):
    return 0.1 * x
