import copy

import numpy as np

_Ax_ext = None
GpuIds = None
def _try_import_Ax_ext():
    global _Ax_ext, GpuIds
    if _Ax_ext is None:
        try:
            from _Ax import _Ax_ext as imported
            _Ax_ext = imported
        except ImportError:
            _Ax_ext = None
    if _Ax_ext is not None and GpuIds is None:
        # GpuIds depends on ctypes as well, so we must import it 
        # after sucessful completion of the _Ax_ext import.
        from .gpu import GpuIds as imported_ids
        GpuIds = imported_ids

def _ensure_Ax_ext_import():
    global _Ax_ext
    if _Ax_ext is None:
        from _Ax import _Ax_ext

_try_import_Ax_ext()


def Ax(img, geo, angles, projection_type="Siddon", **kwargs):
    """Forward projection: simulate the projections of an image.

    Computes :math:`A x` on the GPU, where :math:`A` is the system matrix of
    the scanner described by ``geo`` and ``angles``.

    Parameters
    ----------
    img : numpy.ndarray
        Image, float32, shape ``geo.nVoxel`` (``(z, y, x)``).
    geo : tigre.utilities.geometry.Geometry
        The scanner geometry.
    angles : numpy.ndarray
        Projection angles in radians: shape ``(N,)`` for rotations around z,
        or ``(N, 3)`` for ZYZ Euler angles.
    projection_type : {"Siddon", "interpolated"}, optional
        ``"Siddon"`` (default; ``"ray-voxel"`` is an older name) adds up the
        exact length of each ray inside each voxel. ``"interpolated"`` samples
        the image along each ray every ``geo.accuracy`` voxels, with trilinear
        interpolation.
    gpuids : tigre.utilities.gpu.GpuIds, optional
        GPUs to use. Default: all available.

    Returns
    -------
    numpy.ndarray
        Projections, float32, shape ``(N, *geo.nDetector)`` (``(angle, v, u)``).

    Raises
    ------
    TypeError
        If ``img`` is not float32.
    ValueError
        If ``img`` does not have shape ``geo.nVoxel``, or is complex.

    See Also
    --------
    Atb : backprojection.

    Notes
    -----
    ``Atb`` is not the exact transpose of ``Ax``; see
    :doc:`/user-guide/projectors`.

    Examples
    --------
    >>> import numpy as np
    >>> import tigre
    >>> from tigre.utilities import sample_loader
    >>> geo = tigre.geometry_default(high_resolution=False)
    >>> angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
    >>> img = sample_loader.load_head_phantom(geo.nVoxel)
    >>> proj = tigre.Ax(img, geo, angles)
    >>> proj.shape
    (100, 128, 128)
    """
    _ensure_Ax_ext_import() # check if the import was successful, if not this will raise the original import error with the full message.
    if img.dtype != np.float32:
        raise TypeError("Input data should be float32, not " + str(img.dtype))
    if not np.isreal(img).all():
        raise ValueError("Complex types not compatible for projection.")
    if any(img.shape != geo.nVoxel):
        raise ValueError(
            "Input data should be of shape geo.nVoxel: "
            + str(geo.nVoxel)
            + " not:"
            + str(img.shape)
        )
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    """
    Here we cast all values in geo to single point precision float. This way we know what behavior
    to expect from pytigre to Cuda and can change single parameters accordingly.
    """
    geox.cast_to_single()
    # geox.checknans()

    if kwargs.get("gpuids", None) is None:
        gpuids = GpuIds()
    else:
        gpuids = kwargs["gpuids"]
    # if we have more GPUs than projections to compute, reduce the amount of GPUs. 
    if len(angles) < len(gpuids):
        gpuids.devices = list(gpuids.devices[0:len(angles)])

    return _Ax_ext(img, geox, geox.angles, projection_type, geox.mode, gpuids=gpuids)
