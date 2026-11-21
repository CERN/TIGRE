"""Siddon Ax gives the same result when the image is split in pieces (#859).

When the image does not fit on the GPU, cone-beam Ax projects it in pieces
along z and adds the results. Siddon got some rays wrong by up to ~20%: when a
ray entered a piece through its top or bottom face at a point that was also a
voxel boundary of another axis, the first voxel and the first boundary
crossings disagreed. That is rare for the whole image but common for pieces.

Here the GPU is filled (with plain CUDA memory, no PyTorch needed) so TIGRE
must split a 512^3 image, and the result is compared with the unsplit one. The
pieces add up in a different order, so they agree up to rounding (~1e-5), not
bit for bit.
"""
import ctypes
import os
import subprocess
import sys

import numpy as np
import pytest

import tigre
from tigre.utilities.gpu import GpuIds

pytestmark = [pytest.mark.gpu, pytest.mark.extended]  # extended: not run by a plain `pytest`, see contributing.md


def _cudart():
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
            lib.cudaMalloc.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_size_t]
            lib.cudaFree.argtypes = [ctypes.c_void_p]
            return lib
    return None


CUDART = _cudart()


def free_bytes():
    free, total = ctypes.c_size_t(0), ctypes.c_size_t(0)
    assert CUDART.cudaMemGetInfo(ctypes.byref(free), ctypes.byref(total)) == 0
    return free.value


def fill_gpu(leave, tolerance=8 * 2**20):
    """cudaMalloc so that about `leave` bytes stay free; returns the pointer.

    Measures what is free after allocating and corrects the size: on Windows,
    with other programs drawing on the GPU, the free memory reported before
    allocating can be off by a lot."""
    size = free_bytes() - leave
    for _ in range(4):
        ptr = ctypes.c_void_p()
        if size <= 0 or CUDART.cudaMalloc(ctypes.byref(ptr), size) != 0:
            break
        now = free_bytes()
        if abs(now - leave) <= tolerance:
            return ptr
        CUDART.cudaFree(ptr)
        size += now - leave
    pytest.skip("could not leave {} MB free on the GPU (other programs using it?)".format(leave >> 20))


@pytest.mark.skipif(CUDART is None, reason="CUDA runtime not loadable via ctypes")
def test_siddon_split_matches_unsplit():
    n = 512
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([n, n, n]))
    angles = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    img = np.random.default_rng(0).random((n, n, n), dtype=np.float32)
    one_gpu = GpuIds()
    one_gpu.devices = [0]

    # TIGRE uses 95% of the free memory; minus its 4 blocks of 9 projections,
    # what is left must hold 1/3 to 1/2 of the image: 3 pieces, ~40 MB from
    # either end of that range.
    mem_proj = int(np.prod(geo.nDetector)) * 4
    leave = int((36 * mem_proj + img.nbytes * 5 / 12) / 0.95)
    if free_bytes() < leave + 2**30:
        pytest.skip("not enough free GPU memory for this test")

    whole = tigre.Ax(img, geo, angles, "Siddon", gpuids=one_gpu)
    ptr = fill_gpu(leave)
    try:
        split = tigre.Ax(img, geo, angles, "Siddon", gpuids=one_gpu)
    except RuntimeError as e:
        if "out of memory" in str(e) and free_bytes() < leave - 16 * 2**20:
            pytest.skip("another program took GPU memory during the test")
        raise
    finally:
        CUDART.cudaFree(ptr)

    diff = np.abs(split - whole) / np.abs(whole).max()
    # Before the fix: up to ~0.2. The pieces add up in a different order, and
    # each ray's sum is accurate to ~1e-5 (float32 emulation: ~3e-6).
    assert diff.max() < 1e-4, "max {:.2e} at (angle, v, u) = {}, {:.2%} of pixels above 1e-4".format(
        diff.max(), np.unravel_index(diff.argmax(), diff.shape), (diff > 1e-4).mean())


SYCL_SCRIPT = r"""
import sys
import numpy as np
import tigre

n = 256
geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([n, n, n]))
angles = np.linspace(0, 2 * np.pi, 64, endpoint=False)
img = np.random.default_rng(0).random((n, n, n), dtype=np.float32)
np.save(sys.argv[1], tigre.Ax(img, geo, angles, "Siddon"))
"""


def _sycl_build():
    try:
        from tigre.utilities.cuda_interface import _gpuUtils
    except ImportError:
        return False
    return getattr(_gpuUtils, "getBackend", lambda: "cuda")() == "sycl"


@pytest.mark.skipif(not _sycl_build(), reason="SYCL build only (TIGRE_SYCL_MEMORY_LIMIT_MB)")
def test_siddon_split_matches_unsplit_sycl(tmp_path):
    # SYCL build: TIGRE_SYCL_MEMORY_LIMIT_MB caps the memory TIGRE plans with, so the
    # 64 MB image is projected in pieces on any device. Each run is a subprocess, since
    # the limit is read from the environment.
    script = tmp_path / "siddon_split.py"
    script.write_text(SYCL_SCRIPT)
    out = {}
    for limit in (None, 24):
        env = dict(os.environ)
        env.pop("TIGRE_SYCL_MEMORY_LIMIT_MB", None)
        if limit is not None:
            env["TIGRE_SYCL_MEMORY_LIMIT_MB"] = str(limit)
        path = tmp_path / "ax_{}.npy".format(limit)
        subprocess.run([sys.executable, str(script), str(path)], env=env, check=True, timeout=600)
        out[limit] = np.load(path)
    whole, split = out[None], out[24]
    diff = np.abs(split - whole) / np.abs(whole).max()
    # Before the fix: up to ~0.16 here.
    msg = "max {:.2e} at (angle, v, u) = {}, {:.2%} of pixels above 1e-4"
    assert diff.max() < 1e-4, msg.format(
        diff.max(), np.unravel_index(diff.argmax(), diff.shape), (diff > 1e-4).mean()
    )
