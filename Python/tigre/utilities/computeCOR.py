# -*- coding: utf-8 -*-
"""Centre-of-rotation estimate from a single sinogram slice.

Python port of MATLAB/Utilities/computeCOR.m, which is itself modified from
SophiaBeads
(https://github.com/Sophilyplum/sophiabeads-datasets/blob/master/tools/centre_geom.m)
Reference: T. Liu, "Direct central ray determination in computed
microtomography", Optical Engineering, April 2009.

For a candidate centre of rotation, every ray has a conjugate ray half a turn
away that should measure the same line integral. The candidate that makes the
sinogram agree with its own reflected self is the centre.

Checked against MATLAB/Utilities/computeCOR.m line by line. The arrays are
transposed relative to it throughout - MATLAB holds (detector, angle) and
projections as (V, U, angle), Python holds (angle, detector) and projections
as (angle, V, U) - so the angular wrap is a concatenation along axis 0 here
against a horizontal one there, and `angles[:, None] + beta[None, :]` is the
transpose of its `repmat(angles', nDet, 1) + repmat(beta, 1, nAngles)`.
Everything else - the 21-point search grid per precision, the seven precision
levels, linear interpolation with 0 outside, and `-midpoint * DSO / DSD` -
matches.

One deliberate divergence in the score normalisation, matching MATLAB's
INTENT: MATLAB's `nonzero = find(test > 0)` yields INDICES, so its
`length(nonzero)` is the count of valid pixels and the score is a mean over
them. A boolean-mask translation (`len(mask)` = the angle count) would be
constant across candidates, turning the search into a SUM minimisation that
rewards candidates whose reflected rays fall off the detector - every lost
pixel drops a positive term. The mean over valid pixels is used here.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def computeCOR(data, geo, angles, slc=None):
    """Estimate the centre of rotation, in mm at the rotation axis.

    :param data: projections, (n_angles, n_rows, n_cols)
    :param geo: tigre geometry
    :param angles: projection angles, radians
    :param slc: detector row to use (default: the central row)
    :return: float, the COR offset in mm - the sign convention TIGRE's
        `geo.COR` expects
    """
    if slc is None:
        slc = int(np.floor(data.shape[1] / 2))

    # float64 throughout, as in the MATLAB reference (`double(repmat(...))`):
    # this is one detector row, so the precision is free.
    sino = np.squeeze(np.asarray(data)[:, slc, :]).astype(np.float64)
    angles = np.asarray(angles).astype(np.float64).ravel()

    # check_geo() broadcasts these to one value per angle; the fan-angle
    # geometry below is a single fixed source-detector pair.
    DSD = float(np.mean(np.atleast_1d(np.asarray(geo.DSD, dtype=np.float64))))
    DSO = float(np.mean(np.atleast_1d(np.asarray(geo.DSO, dtype=np.float64))))

    # RegularGridInterpolator needs ascending coordinates; tolerate callers
    # whose angle arrays run backwards.
    if float(angles[0]) > float(angles[-1]):
        angles = angles[::-1]
        sino = sino[::-1]

    det = np.linspace(-geo.sDetector[1] / 2 + geo.dDetector[1] / 2,
                      +geo.sDetector[1] / 2 - geo.dDetector[1] / 2,
                      int(geo.nDetector[1]), dtype=np.float64)
    gamma = np.arctan(det / DSD)                # fan angle per column; candidate-invariant

    # Repeat the sinogram either side in angle so a conjugate ray that lands
    # outside [angles[0], angles[-1]] still interpolates instead of vanishing.
    two_pi = 2.0 * np.pi
    ang3 = np.concatenate((angles - two_pi, angles, angles + two_pi))
    sino3 = np.concatenate((sino, sino, sino), axis=0)
    interp = RegularGridInterpolator((ang3, det), sino3, method="linear",
                                     bounds_error=False, fill_value=0.0)

    midpoint = 0.0
    for pr in (1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001):
        COR = np.linspace(midpoint - 10 * pr, midpoint + 10 * pr, 21, dtype=np.float64)
        M = np.full((COR.size,), np.inf, dtype=np.float64)

        for j in range(COR.size):
            gamma_c = np.arctan(COR[j] / DSD)
            beta = 2.0 * (gamma - gamma_c) + np.pi          # conjugate-ray angle offset
            s2 = DSD * np.tan(2.0 * gamma_c - gamma)        # ... and its detector position

            angles_aux = angles[:, None] + beta[None, :]
            pts = np.stack((angles_aux,
                            np.broadcast_to(s2[None, :], angles_aux.shape)), axis=-1)
            test = interp(pts.reshape(-1, 2)).reshape(angles_aux.shape)

            valid = test > 0
            n_valid = int(valid.sum())
            if n_valid == 0:
                continue        # leaves inf: a candidate that sees nothing never wins
            M[j] = np.sum((test[valid] - sino[valid]) ** 2) / n_valid

        midpoint = float(COR[int(np.argmin(M))])

    return -midpoint * DSO / DSD
