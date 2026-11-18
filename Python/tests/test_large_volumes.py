"""Backprojection of an image with more than 2^31 voxels.

Sizes were computed as products of `int`s in places, which wrap above 2^31:
- matched Atb (`voxel_backprojection2.cu`): the final scaling looped up to
  nVoxelX*nVoxelY*nVoxelZ, an int product;
- parallel-beam Atb (`voxel_backprojection_parallel.cu`): the image buffer was
  allocated with an int product of the three sizes.

A 1300^3 image (2.2e9 voxels, 8.8 GB) is backprojected from 8 projections of
ones on one GPU, unsplit. The geometry is symmetric about the central z plane,
so the result must be too: the first slice equals the last, and so on. Any
part of the image left unscaled or unwritten breaks that.

Needs about 10 GB free on the GPU and 20 GB of RAM; skipped otherwise. Only
when the whole image fits on the GPU does the size reach the code above; with
less GPU memory the image is split in pieces and the test would not show
anything, so it is skipped too.
"""
import ctypes
import os
import sys

import numpy as np
import pytest

import tigre
from tigre.utilities.gpu import GpuIds

pytestmark = pytest.mark.gpu

N = 1300  # 1300^3 = 2.197e9 > 2^31 - 1
N_ANGLES = 8


def _free_gpu_bytes():
    names = ["cudart64_13.dll", "cudart64_12.dll", "cudart64_110.dll",
             "libcudart.so", "libcudart.so.13", "libcudart.so.12", "libcudart.so.11.0"]
    dirs = [""]
    if os.environ.get("CUDA_PATH"):
        dirs.append(os.path.join(os.environ["CUDA_PATH"], "bin"))
    for d in dirs:
        for n in names:
            try:
                lib = ctypes.CDLL(os.path.join(d, n) if d else n)
            except OSError:
                continue
            lib.cudaMemGetInfo.argtypes = [ctypes.POINTER(ctypes.c_size_t)] * 2
            free, total = ctypes.c_size_t(0), ctypes.c_size_t(0)
            if lib.cudaMemGetInfo(ctypes.byref(free), ctypes.byref(total)) == 0:
                return free.value
    return None


def _available_ram_bytes():
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        pass
    if sys.platform == "win32":
        class MemoryStatus(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        status = MemoryStatus()
        status.dwLength = ctypes.sizeof(MemoryStatus)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return status.ullAvailPhys
        return None
    try:
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        return None


def _geometry(mode):
    if mode == "cone":
        return tigre.geometry(mode="cone", default=True, nVoxel=np.array([N, N, N]))
    return tigre.geometry(mode="parallel", nVoxel=np.array([N, N, N]))


@pytest.mark.parametrize("mode,kind", [("cone", "matched"), ("parallel", "FDK")])
def test_backprojection_of_more_than_2_31_voxels(mode, kind):
    image_bytes = N ** 3 * 4
    ram = _available_ram_bytes()
    if ram is None or ram < 2 * image_bytes + 2 * 2**30:
        pytest.skip("needs about {:.0f} GB of free RAM".format((2 * image_bytes + 2 * 2**30) / 2**30))
    free = _free_gpu_bytes()
    if free is None or free < image_bytes + 1.5 * 2**30:
        pytest.skip("needs about {:.0f} GB free on the GPU, so that the image is not split".format(
            (image_bytes + 1.5 * 2**30) / 2**30))

    geo = _geometry(mode)
    angles = np.linspace(0, 2 * np.pi, N_ANGLES, endpoint=False)
    proj = np.ones((N_ANGLES, int(geo.nDetector[0]), int(geo.nDetector[1])), dtype=np.float32)
    one_gpu = GpuIds()
    one_gpu.devices = [0]

    img = tigre.Atb(proj, geo, angles, kind, gpuids=one_gpu)

    assert img.shape == (N, N, N)
    scale = float(np.abs(img[N // 2]).max())
    assert scale > 0
    for k in [0, 1, N // 4, N // 2 - 1]:
        assert np.isfinite(img[k]).all() and np.isfinite(img[N - 1 - k]).all()
        np.testing.assert_allclose(img[k], img[N - 1 - k], rtol=1e-4, atol=1e-5 * scale,
                                   err_msg="slices {} and {} differ".format(k, N - 1 - k))
