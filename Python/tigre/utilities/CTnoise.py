import numpy as np

from .gpu import _import_extension

RNG = _import_extension("_RandomNumberGenerator")

def add(projections, Gaussian=None, Poisson=None, seed=None):
    """Add realistic noise to simulated projections.

    Converts the line integrals to photon counts, adds Poisson (quantum) noise
    and Gaussian (electronic) noise, and converts back. Available as
    ``tigre.utilities.CTnoise.add``.

    Parameters
    ----------
    projections : numpy.ndarray
        Noise-free projections (line integrals).
    Gaussian : numpy.ndarray, optional
        Mean and standard deviation of the electronic noise, in photon counts,
        shape ``(2,)``. Default ``[0, 0.5]``.
    Poisson : float, optional
        Photon count of an unattenuated ray (I0). Lower values give noisier
        data. Default 60000, as MATLAB's ``addCTnoise``.
    seed : int, optional
        Seed of the GPU random number generator: the same seed gives the same
        noise. Default: a new one for every call, drawn from NumPy's global
        random state (so ``np.random.seed(...)`` makes a script reproducible).
        Before TIGRE 3.3 the seed was always the same, so every call added
        exactly the same noise.

    Returns
    -------
    numpy.ndarray
        Noisy projections, float32, same shape.

    Examples
    --------
    >>> from tigre.utilities import CTnoise
    >>> noisy = CTnoise.add(proj, Poisson=1e5, Gaussian=np.array([0, 10]))
    """

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

    projections = RNG.add_noise(projections, Gaussian[0], Gaussian[1], seed=seed)
    # the Gaussian term can push a reading to zero or below, where the log is not finite
    projections[projections <= 0] = 1e-6

    projections = -np.log(projections / Poisson) * max_proj
    projections = np.float32(projections)
    return projections
