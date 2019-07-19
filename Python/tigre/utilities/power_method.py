from __future__ import division
from __future__ import print_function
import numpy as np
import tigre

def svd_power_method(arr,geo,angles,**kwargs):
    maxiter = 100
    if 'maxiter' in kwargs:
        maxiter = kwargs['maxiter']
    #very optimistic value for epsilon (it turns out)
    epsilon = 0.001
    if 'epsilon' in kwargs:
        epsilon = kwargs['epsilon']
    verbose = False
    if 'verbose' in kwargs:
        verbose = kwargs['verbose']

    error = np.inf

    arr_old = arr.copy()
    i = 0
    while error >= epsilon:
        if verbose:
            if i ==0: print('using SVD power method to calculate '
                            'Lipschitz const')
            if np.mod(i,10) ==0:
                print('error for val is:' +str(error))
        # very verbose and inefficient for now
        omega = tigre.Ax(arr_old,geo,angles)
        alpha = np.linalg.norm(omega)
        u = (1./alpha)*omega
        z = tigre.Atb(u,geo,angles)
        beta = np.linalg.norm(z)
        arr = (1./beta)*z
        error = np.linalg.norm(tigre.Ax(arr,geo,angles)-beta*u)
        sigma = beta
        arr_old = arr
        i += 1
        if i >= maxiter:
            return sigma

    return sigma

# http://www.anstuocmath.ro/mathematics/anale2015vol2/Bentbib_A.H.__Kanber_A..pdf