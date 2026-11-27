import copy

import numpy as np
from tigre.utilities.cuda_interface._Atb import _Atb_ext

from . import _gpu_arrays
from .gpu import GpuIds


def Atb(projections, geo, angles, backprojection_type="FDK", **kwargs):
    """Backprojection: spread projections back over the image.

    Approximates :math:`A^T b` on the GPU, where :math:`A` is the system matrix
    of the scanner described by ``geo`` and ``angles``.

    Parameters
    ----------
    projections : numpy.ndarray, torch.Tensor or cupy.ndarray
        Projections, float32, shape ``(N, *geo.nDetector)`` (``(angle, v, u)``).
        A PyTorch CUDA tensor or a CuPy array stays on its GPU: see Notes.
    geo : tigre.utilities.geometry.Geometry
        The scanner geometry.
    angles : numpy.ndarray
        Projection angles in radians: shape ``(N,)`` for rotations around z,
        or ``(N, 3)`` for ZYZ Euler angles.
    backprojection_type : {"FDK", "matched"}, optional
        ``"FDK"`` (default) uses the distance weights of the FDK algorithm.
        ``"matched"`` uses weights that make it close to the transpose of
        :func:`Ax`; the Krylov algorithms use it.
    gpuids : tigre.utilities.gpu.GpuIds, optional
        GPUs to use. Default: all available, or for a GPU array, its GPU.

    Returns
    -------
    numpy.ndarray, torch.Tensor or cupy.ndarray
        Image, float32, shape ``geo.nVoxel`` (``(z, y, x)``), the same kind of
        array as ``projections``, on the same GPU.

    Raises
    ------
    TypeError
        If ``projections`` is not float32.
    ValueError
        If ``projections`` does not have the expected shape, or is complex.

    See Also
    --------
    Ax : forward projection.

    Notes
    -----
    Neither backprojector is the exact transpose of :func:`Ax`; see
    :doc:`/user-guide/projectors`.

    With a PyTorch CUDA tensor or a CuPy array, the data does not go through
    host memory: TIGRE reads ``projections`` on its GPU and returns a new
    array there. This uses that one GPU only. If ``gpuids`` names other GPUs,
    the projections are copied to the host, backprojected as usual and the
    result copied back.
    """
    if not isinstance(projections, np.ndarray) and _gpu_arrays.is_gpu_array(projections):
        return _gpu_arrays.Atb(projections, geo, angles, backprojection_type,
                               kwargs.get("gpuids", None), Atb)
    if projections.dtype != np.float32:
        raise TypeError("Input data should be float32, not " + str(projections.dtype))
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    """
    Here we cast all values in geo to single point precision float. This way we know what behavior
    to expect from pytigre to Cuda and can change single parameters accordingly.
    """
    geox.cast_to_single()
    # geox.checknans()

    if projections.shape != (geox.angles.shape[0], geox.nDetector[0], geox.nDetector[1]):
        raise ValueError(
            "Expected data shape: "
            + str((geox.angles.shape[0], geox.nDetector[0], geox.nDetector[1]))
            + " not compatible with: "
            + str(projections.shape)
        )

    if kwargs.get("gpuids", None) is None:
        gpuids = GpuIds()
    else:
        gpuids = kwargs["gpuids"]

    # if we have more GPUs than slices to compute, reduce the amount of GPUs. 
    if geo.nVoxel[0]< len(gpuids):
        # Trim a copy: the caller's GpuIds (e.g. an algorithm's self.gpuids) must keep all its GPUs (#855).
        gpuids = copy.deepcopy(gpuids)
        gpuids.devices = list(gpuids.devices[0:geo.nVoxel[0]])

    return _Atb_ext(projections, geox, geox.angles, backprojection_type, geox.mode, gpuids=gpuids)
