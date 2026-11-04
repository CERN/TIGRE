import numpy as np


def l2norm(x, chunk=1 << 24):
    """Euclidean norm of a float32 array, accumulated in float64.

    ``np.linalg.norm`` sums the squares in the array's own dtype, so on large
    float32 arrays the result comes back low.

    One pass, no full-size copy. Returns float32 so no caller's dtype changes:
    these scalars multiply float32 volumes, and Ax/Atb reject float64.
    """
    flat = np.ravel(x)
    if flat.dtype != np.float32:
        return np.linalg.norm(flat, 2)
    total = 0.0
    for k in range(0, flat.size, chunk):
        c = flat[k:k + chunk].astype(np.float64)
        total += float(np.dot(c, c))
    return np.float32(np.sqrt(total))


def inner(a, b, chunk=1 << 24):
    """Inner product of two float32 arrays, accumulated in float64.

    Returns float32, like l2norm, so no caller's dtype changes.
    """
    fa, fb = np.ravel(a), np.ravel(b)
    if fa.dtype != np.float32 or fb.dtype != np.float32:
        return np.dot(fa, fb)
    total = 0.0
    for k in range(0, fa.size, chunk):
        total += float(np.dot(fa[k:k + chunk].astype(np.float64),
                              fb[k:k + chunk].astype(np.float64)))
    return np.float32(total)


def im3DNORM(img, normind, varargin=None):
    """
    % IMAGE3DNORM computes the desired image norm
    %   IMAGE3DNORM(IMG,NORMIND) computes the norm of image IMG using the norm
    %   defined in NORMIND
    %
    %   IMG         A 3D image
    %   NORMIND     {non-zero int, inf, -inf, 'fro', 'nuc'}
    %               'TV': TV norm
    %
    %
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % This file is part of the TIGRE Toolbox
    %
    % Copyright (c) 2015, University of Bath and
    %                     CERN-European Organization for Nuclear Research
    %                     All rights reserved.
    %
    % License:            Open Source under BSD.
    %                     See the full license at
    %                     https://github.com/CERN/TIGRE/license.txt
    %
    % Contact:            tigre.toolbox@gmail.com
    % Codes:              https://github.com/CERN/TIGRE/
    % Coded by:           Ander Biguri
    %--------------------------------------------------------------------------
    """
    if normind is [np.inf, -np.inf, "fro", "nuc"]:
        return np.linalg.norm(img.ravel(), normind)
    if type(normind) is int:
        # The L2 case is the one every caller in this toolbox uses, and the one
        # a float32 accumulator cannot compute at these sizes - see l2norm.
        if normind == 2:
            return l2norm(img)
        return np.linalg.norm(img.ravel(), normind)
    if normind == "TV":
        gx, gy, gz = np.diff(img, axis=2), np.diff(img, axis=1), np.diff(img, axis=0)
        g = np.sum(np.sqrt(gx * gx + gy * gy + gz * gz))
        return g
