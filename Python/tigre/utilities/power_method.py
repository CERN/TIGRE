from __future__ import division
from __future__ import print_function

import numpy as np
import tigre
from tigre.utilities.im3Dnorm import l2norm


def svd_power_method(arr, geo, angles, **kwargs):
    """
    method to determine the largest singular value for the rectangular A
    projection matrix.

    :param arr: (np.ndarray, dtype=np.float32)
        random vector
    :param geo: (tigre.geometry())
        geometry of the set up
    :param angles: (np.ndarray, dtype=np.float32)
        angles of projection
    :param kwargs:
        optional arguments
    :keyword maxiter: (int)
        number of iterations to run unless margin of error epsilon has
        been reached. Default: 100
    :keyword epsilon: (float)
        margin of error for the algorithm. Default: 0.001
    :keyword verbose: (bool)
        print stuff. Default: False

    :return: (float)
        the largest singular value of the A matrix.
    """
    maxiter = 100
    if "maxiter" in kwargs:
        maxiter = kwargs["maxiter"]
    # very optimistic value for epsilon (it turns out)
    epsilon = 0.001
    if "epsilon" in kwargs:
        epsilon = kwargs["epsilon"]
    verbose = False
    if "verbose" in kwargs:
        verbose = kwargs["verbose"]

    error = np.inf

    arr_old = arr.copy()
    i = 0
    while error >= epsilon:
        if verbose:
            if i == 0:
                print("using SVD power method to calculate " "Lipschitz const")
            if np.mod(i, 10) == 0:
                print("error for val is:" + str(error))
        # very verbose and inefficient for now
        omega = tigre.Ax(arr_old, geo, angles)
        alpha = l2norm(omega)
        u = (1.0 / alpha) * omega
        # sigma_max(A) needs A^T: the "matched" backprojector is the (near-)
        # adjoint of Ax (inner-product test: ~4e-5 relative asymmetry). The
        # previous call used Atb's default "FDK"-weighted backprojector,
        # which is NOT A^T (~93% inner-product discrepancy) - any Lipschitz/
        # step-size estimate built on it is for the wrong operator.
        z = tigre.Atb(u, geo, angles, backprojection_type="matched")
        beta = l2norm(z)
        arr = (1.0 / beta) * z
        error = l2norm(tigre.Ax(arr, geo, angles) - beta * u)
        sigma = beta
        arr_old = arr
        i += 1
        if i >= maxiter:
            return sigma

    return sigma


# http://www.anstuocmath.ro/mathematics/anale2015vol2/Bentbib_A.H.__Kanber_A..pdf
