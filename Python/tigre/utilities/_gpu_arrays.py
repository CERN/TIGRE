"""GPU arrays in tigre.Ax and tigre.Atb.

tigre.Ax and tigre.Atb also take a PyTorch CUDA tensor or a CuPy array. The
data then stays on the GPU: TIGRE reads the input from the caller's GPU memory
and writes the result into a new array of the same kind on the same GPU,
without going through host memory.

This works on one GPU, the one the array is on. If ``gpuids`` asks for other
GPUs (for example several), the array is copied to the host and the usual
multi-GPU path runs, then the result is copied back: the result is the same,
only slower.

Neither PyTorch nor CuPy is imported here unless the caller passed one of its
arrays, and NumPy input never reaches this module.
"""
import copy


def _kind(x):
    """"torch" or "cupy" for a GPU array, None for anything else."""
    module = type(x).__module__.split(".", 1)[0]
    if module == "torch":
        return "torch" if getattr(x, "is_cuda", False) else None
    if module == "cupy" and hasattr(x, "__cuda_array_interface__"):
        return "cupy"
    return None


def is_gpu_array(x):
    """True for a PyTorch CUDA tensor or a CuPy array."""
    return _kind(x) is not None


# ---- what differs between PyTorch and CuPy
def _prepare(x, kind, name):
    """Check the dtype, make it C-contiguous; return (array, device index)."""
    if kind == "torch":
        import torch

        if x.dtype != torch.float32:
            raise TypeError("{} should be float32, not {}".format(name, x.dtype))
        x = x.detach().contiguous()
        return x, x.device.index
    import cupy

    if x.dtype != cupy.float32:
        raise TypeError("{} should be float32, not {}".format(name, x.dtype))
    x = cupy.ascontiguousarray(x)
    return x, x.device.id


def _empty_like_on(x, kind, shape):
    if kind == "torch":
        import torch

        return torch.empty(tuple(int(s) for s in shape), dtype=torch.float32, device=x.device)
    import cupy

    with x.device:
        return cupy.empty(tuple(int(s) for s in shape), dtype=cupy.float32)


def _wait_for_producer(x, kind):
    """Wait until the work that wrote ``x`` has finished.

    TIGRE runs on its own CUDA streams, so it must not start reading ``x``
    while the stream that wrote it (PyTorch's or CuPy's current stream) is
    still running."""
    if kind == "torch":
        import torch

        torch.cuda.current_stream(x.device).synchronize()
    else:
        with x.device:
            import cupy

            cupy.cuda.get_current_stream().synchronize()


def _pointer(x, kind):
    return x.data_ptr() if kind == "torch" else x.data.ptr


def _to_host(x, kind):
    if kind == "torch":
        return x.detach().cpu().numpy()
    import cupy

    return cupy.asnumpy(x)


def _from_host(result, like, kind):
    if kind == "torch":
        import torch

        return torch.from_numpy(result).to(like.device)
    import cupy

    with like.device:
        return cupy.asarray(result)


# ---- the device to run on
def _gpu_path_ids(gpuids, device):
    """GpuIds for the GPU path, or None if ``gpuids`` asks for other GPUs."""
    if gpuids is not None:
        devices = [int(d) for d in getattr(gpuids, "devices", gpuids)]
        if devices != [device]:
            return None
    from .gpu import GpuIds

    return GpuIds._from_devices([device])


def _checked_geometry(geo, angles):
    geox = copy.deepcopy(geo)
    geox.check_geo(angles)
    geox.cast_to_single()
    return geox


def _keeping_current_device(x, kind):
    """A context that puts the caller's current CUDA device back afterwards.

    TIGRE calls cudaSetDevice and leaves the last GPU it used as the current
    one. PyTorch and CuPy take their current device from CUDA, so without this
    a call on cuda:1 would make later `device="cuda"` tensors go to GPU 1."""
    if kind == "torch":
        import torch

        return torch.cuda.device(x.device)
    return x.device


# ---- the entry points, called from Ax.py and Atb.py
def Ax(img, geo, angles, projection_type, gpuids, host_Ax):
    with _keeping_current_device(img, _kind(img)):
        return _ax(img, geo, angles, projection_type, gpuids, host_Ax)


def Atb(projections, geo, angles, backprojection_type, gpuids, host_Atb):
    with _keeping_current_device(projections, _kind(projections)):
        return _atb(projections, geo, angles, backprojection_type, gpuids, host_Atb)


def _ax(img, geo, angles, projection_type, gpuids, host_Ax):
    kind = _kind(img)
    img, device = _prepare(img, kind, "Input data")
    if tuple(img.shape) != tuple(int(n) for n in geo.nVoxel):
        raise ValueError("Input data should be of shape geo.nVoxel: "
                         + str(geo.nVoxel) + " not:" + str(tuple(img.shape)))
    ids = _gpu_path_ids(gpuids, device)
    if ids is None:
        result = host_Ax(_to_host(img, kind), geo, angles, projection_type, gpuids=gpuids)
        return _from_host(result, img, kind)

    from _Ax import _Ax_ext_gpu

    geox = _checked_geometry(geo, angles)
    out = _empty_like_on(img, kind, (geox.angles.shape[0], geox.nDetector[0], geox.nDetector[1]))
    _wait_for_producer(img, kind)
    _Ax_ext_gpu(_pointer(img, kind), _pointer(out, kind), geox, geox.angles,
                projection_type, geox.mode, gpuids=ids)
    return out


def _atb(projections, geo, angles, backprojection_type, gpuids, host_Atb):
    kind = _kind(projections)
    projections, device = _prepare(projections, kind, "Input data")
    ids = _gpu_path_ids(gpuids, device)
    if ids is None:
        result = host_Atb(_to_host(projections, kind), geo, angles, backprojection_type,
                          gpuids=gpuids)
        return _from_host(result, projections, kind)

    from _Atb import _Atb_ext_gpu

    geox = _checked_geometry(geo, angles)
    expected = (geox.angles.shape[0], int(geox.nDetector[0]), int(geox.nDetector[1]))
    if tuple(projections.shape) != expected:
        raise ValueError("Expected data shape: " + str(expected)
                         + " not compatible with: " + str(tuple(projections.shape)))
    out = _empty_like_on(projections, kind, geox.nVoxel)
    _wait_for_producer(projections, kind)
    _Atb_ext_gpu(_pointer(projections, kind), _pointer(out, kind), geox, geox.angles,
                 backprojection_type, geox.mode, gpuids=ids)
    return out
