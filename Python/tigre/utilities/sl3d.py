"""Three-dimensional Shepp-Logan head phantom
Variations of three-dimensional Shepp-Logan head phantom.

Copyright
=========
3-clause BSD License
Copyright 2021 SADAKANE, Tomoyuki
https://github.com/tsadakane/sl3d

This code is inspired by the following MATLAB codes:
 * Matthias Schabel (2021). 3D Shepp-Logan phantom (https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom), MATLAB Central File Exchange. Retrieved April 29, 2021.
 * Patrick Bolan (2021). 3D Shepp-Logan Phantom (https://www.mathworks.com/matlabcentral/fileexchange/50974-3d-shepp-logan-phantom), MATLAB Central File Exchange. Retrieved April 29, 2021.
"""

import numpy as np

#%%
def shepp_logan_3d(size_out=128, phantom_type="yu-ye-wang", get_ellipsoids=False):
    """Three-dimensional Shepp-Logan head phantom
    Variations of three-dimensional Shepp-Logan head phantom.
    Parameters
    ==========
    phantom_type: str
        One of {"kak-slaney", "yu-ye-wang", "toft-schabel"}, optional
        Default is "yu-ye-wang"
        The type of phantom.

    size_out : int or list whose length is three, optional
        Default is [128, 128, 128]
        The number of voxels of phantom. [nVoxelZ, nVoxelY, nVoxelX]

    get_ellipsoids: bool
        Default is False.
        If False, returns the parameters of the ellipsoids as well as the phantom image.

    Returns
    =======
    If get_ellipsoids is True:
        (img_phantom, ellipsoids)
        img_phantom: numpy ndarray
        ellipsoids: list of parameters that defines ellipsoids
    else:
        img_phantom: numpy ndarray

    Notes
    =====
    For any given voxel in the output image, the voxel's value is equal to the
    sum of the additive intensity values of all ellipsoids that the voxel is a
    part of.  If a voxel is not part of any ellipsoid, its value is 0.
    The additive intensity value A for an ellipsoid can be positive or
    negative;  if it is negative, the ellipsoid will be darker than the
    surrounding voxels.
    Note that, depending on the values of A, some voxels may have values
    outside the range [0,1].
    The voxel value at (x, y, z) is img[z, y, x].
    """

    (ellipsoids, nVoxel, formula) = _parse_inputs(size_out, phantom_type)
    nVoxelX = nVoxel[2]
    nVoxelY = nVoxel[1]
    nVoxelZ = nVoxel[0]

    img_phantom = np.zeros(nVoxel, dtype=np.float32)
    range_x = np.linspace(-1, +1, nVoxelX)
    range_y = np.linspace(-1, +1, nVoxelY)
    range_z = np.linspace(-1, +1, nVoxelZ)
    mesh_z, mesh_y, mesh_x = np.meshgrid(range_z, range_y, range_x, indexing="ij")

    mesh_x = mesh_x.reshape(-1)
    mesh_y = mesh_y.reshape(-1)
    mesh_z = mesh_z.reshape(-1)

    coord = np.vstack([mesh_z, mesh_y, mesh_x])
    img_phantom = img_phantom.reshape(-1)
    for ellipsoid in ellipsoids:
        asq = ellipsoid[0] ** 2  # a^2
        bsq = ellipsoid[1] ** 2  # b^2
        csq = ellipsoid[2] ** 2  # c^2
        x0 = ellipsoid[3]  # x offset
        y0 = ellipsoid[4]  # y offset
        z0 = ellipsoid[5]  # z offset
        phi1 = ellipsoid[6] * np.pi / 180  # 1st Euler angle in radians (rotation about z-axis)
        phi2 = ellipsoid[7] * np.pi / 180  # 2nd Euler angle in radians (rotation about x'-axis)
        phi3 = ellipsoid[8] * np.pi / 180  # 3rd Euler angle in radians (rotation about z"-axis)
        A = ellipsoid[9]  # Amplitude change for this ellipsoid

        c1 = np.cos(phi1)
        s1 = np.sin(phi1)
        c2 = np.cos(phi2)
        s2 = np.sin(phi2)
        c3 = np.cos(phi3)
        s3 = np.sin(phi3)
        # Euler rotation matrix
        alpha = [  # Z      Y                    X
            [c2, -s2 * c1, s2 * s1],  # Z
            [c3 * s2, -s3 * s1 + c3 * c2 * c1, -s3 * c1 - c3 * c2 * s1],  # Y
            [s3 * s2, c3 * s1 + s3 * c2 * c1, c3 * c1 - s3 * c2 * s1],  # X
        ]
        if formula == 0:
            # Move the ellipsoid to the origin first, and rotate...
            coord_rot = np.dot(alpha, coord - np.array([[z0], [y0], [x0]]))
            idx = np.argwhere(
                (coord_rot[2, :]) ** 2 / asq
                + (coord_rot[1, :]) ** 2 / bsq
                + (coord_rot[0, :]) ** 2 / csq
                <= 1
            )
            # Naive:
            # coord_rot = np.dot(alpha, coord-np.array([[z0], [y0], [x0]])) + np.array([[z0], [y0], [x0]])
            # idx = np.argwhere((coord_rot[2,:]-x0)**2/asq + (coord_rot[1,:]-y0)**2/bsq + (coord_rot[0,:]-z0)**2/csq <= 1)
        else:
            # (x0,y0,z0) rotates too!
            coord_rot = np.dot(alpha, coord)
            idx = np.argwhere(
                (coord_rot[2, :] - x0) ** 2 / asq
                + (coord_rot[1, :] - y0) ** 2 / bsq
                + (coord_rot[0, :] - z0) ** 2 / csq
                <= 1
            )
        img_phantom[idx] += A

    img_phantom = img_phantom.reshape(nVoxel)
    if get_ellipsoids:
        return img_phantom, ellipsoids
    else:
        return img_phantom


def _parse_inputs(size_out, phantom_type):
    """
    Returns:
     tuple (ellipsoids, nVoxel)
     * ellipsoids is the m-by-10 array which defines m ellipsoids,
       where m is 10 in the cases of the variants implemented in this file.
     * nVoxel is the 3 array which defines the number of voxels
    Parameters:
     phantom_type: One of {"kak-slaney", "yu-ye-wang", "toft-schabel"}
     size_out: An int or 3-vector.
       * int : the phantom voxel will be isotropic.
       * 3-vector: the size of the phantom image [nZ, nY, nX]
    """
    if type(size_out) == int:
        nVoxel = [size_out, size_out, size_out]
    elif (type(size_out) == list or type(size_out) == tuple) and len(size_out) == 3:
        nVoxel = [size_out[0], size_out[1], size_out[2]]
    elif type(size_out) == np.array and np.size(size_out) == 3:
        nVoxel = [size_out.reshape(-1)[0], size_out.reshape(-1)[1], size_out.reshape(-1)[2]]
    else:
        nVoxel = [128, 128, 128]

    if phantom_type == "kak-slaney":
        ellipsoids = kak_slaney()
        formula = 0
    elif phantom_type == "yu-ye-wang":
        ellipsoids = yu_ye_wang()
        formula = 0
    elif phantom_type == "toft-schabel":
        ellipsoids = toft_schabel()
        formula = 1
    else:
        print(f"Unknown type {phantom_type}. yu-ye-wang is used.")
        ellipsoids = yu_ye_wang()
        formula = 0

    return (ellipsoids, nVoxel, formula)


###################################
#  Definetions of Head phantoms:  #
###################################
def kak_slaney():
    """
    The 3D Shepp-Logan head phantom. A is the relative density of water.
    Ref:
     [1] Kak AC, Slaney M, Principles of Computerized Tomographic Imaging, 1988. p.102
         http://www.slaney.org/pct/pct-errata.html
    """
    #            a        b        c      x0       y0      z0    phi1  phi2   phi3   A
    #        -------------------------------------------------------------------------
    ells = [
        [0.6900, 0.920, 0.900, 0.000, 0.000, 0.000, 0, 0, 0, 2.00],
        [0.6624, 0.874, 0.880, 0.000, 0.000, 0.000, 0, 0, 0, -0.98],
        [0.4100, 0.160, 0.210, -0.220, 0.000, -0.250, 108, 0, 0, -0.02],
        [0.3100, 0.110, 0.220, 0.220, 0.000, -0.250, 72, 0, 0, -0.02],
        [0.2100, 0.250, 0.500, 0.000, 0.350, -0.250, 0, 0, 0, 0.02],
        [0.0460, 0.046, 0.046, 0.000, 0.100, -0.250, 0, 0, 0, 0.02],
        [0.0460, 0.023, 0.020, -0.080, -0.650, -0.250, 0, 0, 0, 0.01],
        [0.0460, 0.023, 0.020, 0.060, -0.650, -0.250, 90, 0, 0, 0.01],
        [0.0560, 0.040, 0.100, 0.060, -0.105, 0.625, 90, 0, 0, 0.02],
        [0.0560, 0.056, 0.100, 0.000, 0.100, 0.625, 0, 0, 0, -0.02],
    ]
    ells = np.asarray(np.matrix(ells))
    return ells


def yu_ye_wang():
    """
    A variant of the Kak-Slaney phantom in which the contrast is improved for better
    visual perception.
    Ref:
     [2] Yu H, Ye Y, Wang G, Katsevich-Type Algorithms for Variable Radius Spiral Cone-Beam CT
         Proceedings of the SPIE, Volume 5535, p. 550-557 (2004)
    """
    ells = kak_slaney()
    ells[:, 9] = np.array([1.00, -0.80, -0.2, -0.2, 0.2, 0.2, 0.1, 0.1, 0.2, -0.2])
    # #            a      b       c      x0       y0      z0   phi1 phi2 phi3   A
    # #        -----------------------------------------------------------------
    # ells = [[ 0.6900,  0.920,  0.900,   0   ,   0    ,  0    , 0  , 0, 0,   1.0 ],
    #         [ 0.6624,  0.874,  0.880,   0   ,   0    ,  0    , 0  , 0, 0,  -0.8 ],
    #         [ 0.4100,  0.160,  0.210,  -0.22,   0    , -0.250, 108, 0, 0,  -0.2 ],
    #         [ 0.3100,  0.110,  0.220,   0.22,   0    , -0.25 , 72 , 0, 0,  -0.2 ],
    #         [ 0.2100,  0.250,  0.500,   0   ,   0.35 , -0.25 , 0  , 0, 0,   0.2 ],
    #         [ 0.0460,  0.046,  0.046,   0   ,   0.1  , -0.25 , 0  , 0, 0,   0.2 ],
    #         [ 0.0460,  0.023,  0.020,  -0.08,  -0.65 , -0.25 , 0  , 0, 0,   0.1 ],
    #         [ 0.0460,  0.023,  0.020,   0.06,  -0.65 , -0.25 , 90 , 0, 0,   0.1 ],
    #         [ 0.0560,  0.040,  0.100,   0.06,  -0.105,  0.625, 90 , 0, 0,   0.2 ],
    #         [ 0.0560,  0.056,  0.100,   0   ,   0.100,  0.625, 0  , 0, 0,  -0.2 ]]
    return ells


def toft_schabel():
    """
    The geometry of this phantom is based on the 2D phantom shown in [3] and [4].
    (Maybe this is the original Shepp-Logan head phantom?)
    In [5], the intensities of the Shepp-Logan are modified
    to yield higher contrast in the image.
    It is known as 'Modified Shepp-Logan' of the `phantom` function of "Image Processing Toolbox" for MATLAB
    In [6], it is extended to the 3D version. The parameters are as below.
    The formula of geometry transfom for this option is the same as of [6] to reproduce the result,
    while for other options, kak-slaney and yu-ye-wang, it is different.

    Ref:
     [3] Kak AC, Slaney M, Principles of Computerized Tomographic Imaging, 1988. p.55
     [4] Jain, Anil K., Fundamentals of Digital Image Processing, Englewood Cliffs, NJ, Prentice Hall, 1989, p. 439.
     [5] Toft P, The Radon Transform: Theory and Implementation, 1996.
     [6] Matthias Schabel (2021). 3D Shepp-Logan phantom
       (https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom),
       MATLAB Central File Exchange. Retrieved April 29, 2021.
    """
    #              a     b     c     x0      y0      z0    phi1  phi2   phi3   A
    #        -----------------------------------------------------------------
    ells = [
        [0.6900, 0.9200, 0.810, 0, 0, 0, 0, 0, 0, 1.0],
        [0.6624, 0.8740, 0.780, 0, -0.0184, 0, 0, 0, 0, -0.8],
        [0.1100, 0.3100, 0.220, 0.22, 0, 0, -18, 0, 10, -0.2],
        [0.1600, 0.4100, 0.280, -0.22, 0, 0, 18, 0, 10, -0.2],
        [0.2100, 0.2500, 0.410, 0, 0.35, -0.15, 0, 0, 0, 0.1],
        [0.0460, 0.0460, 0.050, 0, 0.1, 0.25, 0, 0, 0, 0.1],
        [0.0460, 0.0460, 0.050, 0, -0.1, 0.25, 0, 0, 0, 0.1],
        [0.0460, 0.0230, 0.050, -0.08, -0.605, 0, 0, 0, 0, 0.1],
        [0.0230, 0.0230, 0.020, 0, -0.606, 0, 0, 0, 0, 0.1],
        [0.0230, 0.0460, 0.020, 0.06, -0.605, 0, 0, 0, 0, 0.1],
    ]
    return ells


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    nVoxelZ = 64  # 256
    nVoxelY = 256  # 128
    nVoxelX = 256  # 64

    phantom0, e = shepp_logan_3d(
        size_out=[nVoxelZ, nVoxelY, nVoxelX], phantom_type="kak-slaney", get_ellipsoids=True
    )
    phantom1 = shepp_logan_3d(size_out=[nVoxelZ, nVoxelY, nVoxelX], phantom_type="yu-ye-wang")
    phantom2 = shepp_logan_3d(size_out=[nVoxelZ, nVoxelY, nVoxelX], phantom_type="toft-schabel")
    plt.imshow(
        np.concatenate(
            [
                phantom0[3 * phantom0.shape[0] // 8],
                phantom1[3 * phantom1.shape[0] // 8],
                phantom2[phantom2.shape[0] // 2],
            ],
            axis=1,
        )
    )
    plt.colorbar()
    plt.title("KS z=-0.25       YYW z=-0.25        TS z=0.00")
    plt.show()
    print(f"Ellipsoids of kak-slaney are\n{e}")
