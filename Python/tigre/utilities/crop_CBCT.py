import numpy as np

def cropCBCT(img):
    """
    Apply cylindrical CBCT FOV mask to a reconstructed volume.

    Parameters
    ----------
    img : ndarray
        Reconstructed volume in (Z, Y, X) order.

    Returns
    -------
    img_masked : ndarray
        Volume with voxels outside CBCT FOV set to zero.
    """
    Nz, Ny, Nx = img.shape
    cx, cy = Nx // 2, Ny // 2
    radius = min(cx, cy)

    Y, X = np.ogrid[:Ny, :Nx]
    dist = np.sqrt((X - cx)**2 + (Y - cy)**2)
    mask2d = dist <= radius

    img_masked = img.copy()
    img_masked[:, ~mask2d] = 0.0

    return img_masked
