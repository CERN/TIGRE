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
    %   NORMIND     Norm order for the flattened image vector:
    %               - 2: Euclidean norm (L2, optimized for float32 arrays)
    %               - int or float (e.g. 1, 0, -1, -2, 1.5): vector p-norm
    %               - inf, -inf (np.inf, -np.inf): max / min absolute value
    %               - 'TV': Total Variation norm
    %
    %               Note: Matrix norms 'fro' (Frobenius) and 'nuc' (nuclear)
    %               are unsupported because im3DNORM flattens the input array
    %               into a 1D vector (img.ravel()). Passing 'fro' or 'nuc'
    %               will raise a ValueError.
    %
    %   VARARGIN    Optional arguments for specific norms:
    %               - For 'TV': gradient type {'forward' (default), 'backward', 'central'}
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
    if normind is None or isinstance(normind, bool) or not np.isscalar(normind):
        raise ValueError(f"Unsupported norm order: {normind}")
    if isinstance(normind, (float, np.floating)) and np.isnan(normind):
        raise ValueError(f"Unsupported norm order: {normind}")
    if normind == "TV":
        img = np.asarray(img)
        if img.ndim != 3:
            raise ValueError(f"TV norm requires a 3D image, got {img.ndim}D array")
        type_grad = "forward"
        if varargin is not None:
            if isinstance(varargin, (list, tuple)) and len(varargin) > 0:
                type_grad = varargin[0]
            elif isinstance(varargin, str):
                type_grad = varargin
            else:
                raise ValueError(f"Unknown gradient option for TV norm: {varargin}")
        if not isinstance(type_grad, str):
            raise ValueError(f"Unknown gradient option for TV norm: {type_grad}")
        type_grad = type_grad.lower()
        if type_grad == "forward":
            gx = np.pad(np.diff(img, axis=2), ((0, 0), (0, 0), (0, 1)), mode="constant")
            gy = np.pad(np.diff(img, axis=1), ((0, 0), (0, 1), (0, 0)), mode="constant")
            gz = np.pad(np.diff(img, axis=0), ((0, 1), (0, 0), (0, 0)), mode="constant")
        elif type_grad == "backward":
            gx = np.pad(np.diff(img, axis=2), ((0, 0), (0, 0), (1, 0)), mode="constant")
            gy = np.pad(np.diff(img, axis=1), ((0, 0), (1, 0), (0, 0)), mode="constant")
            gz = np.pad(np.diff(img, axis=0), ((1, 0), (0, 0), (0, 0)), mode="constant")
        elif type_grad == "central":
            gz, gy, gx = np.gradient(img)
        else:
            raise ValueError(f"Unknown gradient option for TV norm: {type_grad}")
        g = np.sum(np.sqrt(gx * gx + gy * gy + gz * gz))
        return g
    # The L2 case is the one every caller in this toolbox uses, and the one
    # a float32 accumulator cannot compute at these sizes - see l2norm.
    if normind == 2:
        return l2norm(img)
    try:
        return np.linalg.norm(np.ravel(img), normind)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Unsupported norm order: {normind}") from exc
