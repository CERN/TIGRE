from __future__ import division
from __future__ import print_function
import os
import sys
# TODO: this is quite nasty; it would be nice to reorganise file structure later so top level folder is always in path
currDir = os.path.dirname(os.path.realpath(__file__))
rootDir = os.path.abspath(os.path.join(currDir, '..'))
if rootDir not in sys.path:  # add parent dir to paths
    sys.path.append(rootDir)

import numpy as np

from _Ax import Ax
from _Atb import Atb


def CGLS(projections, geometry, angles, iterations, verbose=False, log_parameters=False):
    if log_parameters:
        parameter_history = {}
        parameter_history['alpha'] = np.zeros([iterations], dtype=np.float32)
        parameter_history['beta'] = np.zeros([iterations], dtype=np.float32)
        parameter_history['gamma'] = np.zeros([iterations], dtype=np.float32)
        parameter_history['q_norm'] = np.zeros([iterations], dtype=np.float32)
        parameter_history['s_norm'] = np.zeros([iterations], dtype=np.float32)

    x = np.zeros(geometry.nVoxel.astype(int), dtype=np.float32)

    r = projections - Ax(x, geometry, angles, 'ray-voxel')
    p = Atb(r, geometry, angles, 'matched')
    p_norm = np.linalg.norm(p.ravel(), 2)
    gamma = p_norm * p_norm

    errorsL2 = np.zeros([iterations], dtype=np.float32)
    for i in range(iterations):
        q = Ax(p, geometry, angles, 'ray-voxel')
        q_norm = np.linalg.norm(q.ravel(), 2)
        alpha = gamma/(q_norm * q_norm)
        x += alpha * p

        aux = projections - Ax(x, geometry, angles, 'ray-voxel')
        errorsL2[i] = np.linalg.norm(aux.ravel(), 2)

        if i > 0 and errorsL2[i] > errorsL2[i-1]:
            x -= alpha * p
            if verbose:
                print('CGLS stoped in iteration', i, 'due to divergence.')
            if log_parameters:
                return x, errorsL2, parameter_history
            return x, errorsL2

        r -= alpha * q

        s = Atb(r, geometry, angles, 'matched')
        s_norm = np.linalg.norm(s.ravel(), 2)

        gamma1 = s_norm * s_norm
        beta = gamma1/gamma

        if log_parameters:
            parameter_history['alpha'][i] = alpha
            parameter_history['beta'][i] = beta
            parameter_history['gamma'][i] = gamma
            parameter_history['q_norm'][i] = q_norm
            parameter_history['s_norm'][i] = s_norm

        gamma = gamma1
        p = s + beta * p

    if log_parameters:
        return x, errorsL2, parameter_history
    return x, errorsL2
