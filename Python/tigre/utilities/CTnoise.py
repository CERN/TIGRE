import numpy as np
import _RandomNumberGenerator as RNG

def add(projections, Gaussian=None, Poisson=None):

    if Poisson is not None:
        if not np.isscalar(Poisson):
            raise ValueError(
                "Poisson value should be an scalar, is " + str(type(Poisson)) + " instead."
            )
    else:
        # same default photon count as MATLAB's addCTnoise
        Poisson = 60000
        if np.max(projections) > Poisson:
            Poisson = np.max(projections) / 5
    if Gaussian is not None:
        if not isinstance(Gaussian, np.ndarray):
            raise ValueError(
                "Gaussian value should be an array, is " + str(type(Gaussian)) + " instead."
            )
        if Gaussian.shape != (2,):
            raise ValueError("Gaussian shape should be 1x2, is " + str(Gaussian.shape) + "instead.")
    else:
        Gaussian = np.array([0, 0.5])
    max_proj = np.max(projections)
    projections = Poisson * np.exp(-projections / max_proj)

    projections = RNG.add_noise(projections, Gaussian[0], Gaussian[1])
    # the Gaussian term can push a reading to zero or below, where the log is not finite
    projections[projections <= 0] = 1e-6

    projections = -np.log(projections / Poisson) * max_proj
    projections = np.float32(projections)
    return projections
