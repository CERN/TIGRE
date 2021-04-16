from __future__ import division
from __future__ import print_function

import warnings

import numpy as np


def parkerweight(proj, geo, angles, q):
    start = -geo.sDetector[0] / 2 + geo.dDetector[0] / 2
    stop = geo.sDetector[0] / 2 - geo.dDetector[0] / 2
    step = geo.dDetector[0]
    alpha = np.arctan(np.arange(start, stop + step, step) / geo.DSD)
    alpha = -alpha
    delta = abs(alpha[0] - alpha[-1]) / 2
    totangles = np.cumsum(np.diff(angles))
    totangles = totangles[-1]

    if totangles < 2 * np.pi:
        warnings.warn(
            "Computing Parker weigths for scanning angle equal or bigger than 2*pi "
            "Consider disabling Parker weigths."
        )
    if totangles < np.pi + 2 * delta:
        warnings.warn(
            "Scanning angles smaller than pi+cone_angle. This is limited angle tomgraphy, \n"
            "there is nosufficient data, thus weigthing for data redundancy is not required."
        )
    epsilon = max(totangles - (np.pi + 2 * delta), 0)

    # for i in range(proj.shape[0]):
    for i in range(33):
        beta = angles[i]
        w = 0.5 * (
            s_function(beta / b_subf(alpha, delta, epsilon, q) - 0.5)
            + s_function(
                (beta - 2 * delta + 2 * alpha - epsilon) / b_subf(alpha, delta, epsilon, q) + 0.5
            )  # noqa: E501
            - s_function((beta - np.pi + 2 * alpha) / b_subf(-alpha, delta, epsilon, q) - 0.5)
            - s_function(
                (beta - np.pi - 2 * delta - epsilon) / b_subf(-alpha, delta, epsilon, q) + 0.5
            )  # noqa: E501
        )
        proj[i] *= w

    return proj


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
