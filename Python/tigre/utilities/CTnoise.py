import numpy as np
import _RandomNumberGenerator as RNG

def add(projections, Gaussian=None, Poisson=None):

    if Poisson is not None:
        if not np.isscalar(Poisson):
            raise ValueError(
                "Poisson value should be an scalar, is " + str(type(Poisson)) + " instead."
            )
    else:
        Poisson = np.ceil(np.log2(np.max(np.abs(projections))))  # nextpow2
    if Gaussian is not None:
        if not isinstance(Gaussian, np.ndarray):
            raise ValueError(
                "Gaussian value should be an array, is " + str(type(Gaussian)) + " instead."
            )
        if Gaussian.shape != (2,):
            raise ValueError("Gaussian shape should be 1x2, is " + str(Gaussian.shape) + "instead.")

    max_proj = np.max(projections)
    projections = Poisson * np.exp(-projections / max_proj)

    projections = RNG.add_noise(projections, Gaussian[0], Gaussian[1])

    projections = -np.log(projections / Poisson) * max_proj
    projections = np.float32(projections)
    return projections
