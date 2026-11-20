from __future__ import division
from __future__ import print_function

import warnings

import numpy as np

# Wesarg, Stefan, Matthias Ebert, and Thomas Bortfeld. 
# "Parker weights revisited." Medical physics 29.3 (2002): 372-378.

# @article{wesarg2002parker,
#   title={Parker weights revisited},
#   author={Wesarg, Stefan and Ebert, Matthias and Bortfeld, Thomas},
#   journal={Medical physics},
#   volume={29},
#   number={3},
#   pages={372--378},
#   year={2002},
#   publisher={American Association of Physicists in Medicine}
# }


def parkerweight(proj, geo, angles, q):
    """Multiply projections by Parker short-scan redundancy weights, in place.

    proj has shape (nangles, V, U); the fan of a circular scan runs along U.
    q is the smoothing parameter of Wesarg et al. (q=1 gives Parker's weights).
    Follows MATLAB's ParkerWeight.m.
    """
    angles = np.asarray(angles, dtype=float)
    if angles.ndim > 1:
        angles = angles[:, 0]
    if np.unique(np.atleast_1d(geo.DSD)).size > 1:
        warnings.warn("Parker weights not supported for varying geo.DSD, mean(geo.DSD) is used.")
    start = -geo.sDetector[1] / 2 + geo.dDetector[1] / 2
    step = geo.dDetector[1]
    us = start + step * np.arange(int(geo.nDetector[1]))
    alpha = np.arctan(us / np.mean(geo.DSD))
    if np.all(np.diff(angles) >= 0):
        # for increasing angles, flip alpha left-right (as MATLAB does)
        alpha = -alpha
    delta = abs(alpha[0] - alpha[-1]) / 2
    totangles = abs(angles[0] - angles[-1])

    if totangles >= 2 * np.pi:
        warnings.warn(
            "Computing Parker weights for scanning angle equal or bigger than 2*pi "
            "Consider disabling Parker weights."
        )
    if totangles < np.pi + 2 * delta:
        warnings.warn(
            "Scanning angles smaller than pi+cone_angle. This is limited angle tomography, \n"
            "there is no sufficient data, thus weighting for data redundancy is not required."
        )
    epsilon = max(totangles - (np.pi + 2 * delta), 0)

    for ii in range(proj.shape[0]):
        beta = abs(angles[ii] - angles[0])
        # at the edge column of a scan of exactly pi + fan angle, b is 0: the 0/0
        # gives NaN, which s_function maps to 0, as in MATLAB
        with np.errstate(divide="ignore", invalid="ignore"):
            w = _weights(beta, alpha, delta, epsilon, q)
        # w varies along U, the last axis
        proj[ii] *= w[np.newaxis, :].astype(proj.dtype)

    return proj


def _weights(beta, alpha, delta, epsilon, q):
    return 0.5 * (
            s_function(beta / b_subf(alpha, delta, epsilon, q) - 0.5)
            + s_function(
                (beta - 2 * delta + 2 * alpha - epsilon) / b_subf(alpha, delta, epsilon, q) + 0.5
            )
            - s_function((beta - np.pi + 2 * alpha) / b_subf(-alpha, delta, epsilon, q) - 0.5)
            - s_function(
                (beta - np.pi - 2 * delta - epsilon) / b_subf(-alpha, delta, epsilon, q) + 0.5
            )
        )


def s_function(abeta):
    w = np.zeros(abeta.shape)
    w[abeta <= -0.5] = 0
    w[abs(abeta) < 0.5] = 0.5 * (1 + np.sin(np.pi * abeta[abs(abeta) < 0.5]))
    w[abeta >= 0.5] = 1
    return w


# b_function = B
def b_function(alpha, delta, epsilon):
    return 2 * delta - 2 * alpha + epsilon


# b_subf = b
def b_subf(alpha, delta, epsilon, q):
    return q * b_function(alpha, delta, epsilon)
