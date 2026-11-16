"""CUDA errors raise a Python exception and leave TIGRE usable (#853).

Before, a CUDA error inside TIGRE called exit(1): the whole Python process
ended, with no exception to catch. Now it raises RuntimeError, and the failed
call releases what it held (device memory, streams, textures, and the host
registration of the caller's array), so the process and TIGRE carry on.

A reliable, cheap CUDA error: a detector (or image) 131072 elements wide. The
texture that holds it would be wider than any GPU supports (at most 16384), so
creating it fails, with almost no memory involved.

Also checks that "GPU heavily used by another program" is a warning, shown
once, instead of an error.
"""
import ctypes
import os
import warnings

import numpy as np
import pytest

import tigre

pytestmark = pytest.mark.gpu

TOO_WIDE = 131072


# ---- direct access to the CUDA runtime, for memory readings and registration checks
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
            lib.cudaHostRegister.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_uint]
            lib.cudaHostUnregister.argtypes = [ctypes.c_void_p]
            return lib
    return None


CUDART = _cudart()
needs_cudart = pytest.mark.skipif(CUDART is None, reason="CUDA runtime not loadable via ctypes")


def mem_info():
    free, total = ctypes.c_size_t(0), ctypes.c_size_t(0)
    assert CUDART.cudaMemGetInfo(ctypes.byref(free), ctypes.byref(total)) == 0
    return free.value, total.value


# ---- problems that fail, and ones that work
def too_wide_detector(mode, nvoxel=8):
    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array([nvoxel] * 3))
    geo.nDetector = np.array([1, TOO_WIDE])
    geo.dDetector = np.array([0.8, 0.8])
    geo.sDetector = geo.nDetector * geo.dDetector
    angles = np.linspace(0, 2 * np.pi, 4, endpoint=False, dtype=np.float32)
    proj = np.ones((len(angles), 1, TOO_WIDE), dtype=np.float32)
    return geo, angles, proj


def too_wide_image(mode):
    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array([1, 1, TOO_WIDE]))
    angles = np.linspace(0, 2 * np.pi, 4, endpoint=False, dtype=np.float32)
    img = np.ones(tuple(geo.nVoxel), dtype=np.float32)
    return geo, angles, img


def small_problem():
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([32, 32, 32]))
    angles = np.linspace(0, 2 * np.pi, 40, endpoint=False, dtype=np.float32)
    img = np.random.default_rng(0).random(tuple(geo.nVoxel), dtype=np.float32)
    return geo, angles, img


def call_failing(kind):
    if kind.startswith("Atb"):
        mode, bp = {"Atb cone FDK": ("cone", "FDK"), "Atb cone matched": ("cone", "matched"),
                    "Atb parallel": ("parallel", "FDK")}[kind]
        geo, angles, proj = too_wide_detector(mode)
        tigre.Atb(proj, geo, angles, bp)
    else:
        mode, pt = {"Ax cone Siddon": ("cone", "Siddon"), "Ax cone interpolated": ("cone", "interpolated"),
                    "Ax parallel Siddon": ("parallel", "Siddon"),
                    "Ax parallel interpolated": ("parallel", "interpolated")}[kind]
        geo, angles, img = too_wide_image(mode)
        tigre.Ax(img, geo, angles, pt)


KINDS = ["Atb cone FDK", "Atb cone matched", "Atb parallel",
         "Ax cone Siddon", "Ax cone interpolated", "Ax parallel Siddon", "Ax parallel interpolated"]


@pytest.mark.parametrize("kind", KINDS)
def test_cuda_error_raises_and_tigre_keeps_working(kind):
    geo, angles, img = small_problem()
    expected_proj = tigre.Ax(img, geo, angles)
    expected_img = tigre.Atb(expected_proj, geo, angles)

    # Reaching the assertion at all is half the test: before, the process exited.
    with pytest.raises(RuntimeError):
        call_failing(kind)

    np.testing.assert_array_equal(tigre.Ax(img, geo, angles), expected_proj)
    np.testing.assert_array_equal(tigre.Atb(expected_proj, geo, angles), expected_img)


@needs_cudart
@pytest.mark.parametrize("mode,bp", [("cone", "FDK"), ("cone", "matched"), ("parallel", "FDK")])
def test_failed_call_releases_the_callers_array(mode, bp):
    """Atb page-locks (registers) the caller's projections. If a failed call
    left them registered, every later call with the same array would fail."""
    geo, angles, proj = too_wide_detector(mode)
    with pytest.raises(RuntimeError):
        tigre.Atb(proj, geo, angles, bp)
    err = CUDART.cudaHostRegister(proj.ctypes.data, proj.nbytes, 0)
    assert err == 0, "projections still registered after the failed call (cudaHostRegister error {})".format(err)
    CUDART.cudaHostUnregister(proj.ctypes.data)


@needs_cudart
def test_failed_calls_release_device_memory():
    """Cone Atb allocates the image on the GPU before it fails on the texture.
    Without the release, each failed call would keep that memory."""
    geo, angles, proj = too_wide_detector("cone", nvoxel=256)  # 64 MB image on the GPU
    with pytest.raises(RuntimeError):
        tigre.Atb(proj, geo, angles)  # absorb one-off setup
    before, _ = mem_info()
    for _ in range(20):
        with pytest.raises(RuntimeError):
            tigre.Atb(proj, geo, angles)
    after, _ = mem_info()
    drift_mb = (before - after) / 2**20
    assert drift_mb < 100, "{:.0f} MB of GPU memory lost over 20 failed calls".format(drift_mb)


@needs_cudart
def test_busy_gpu_is_a_warning_shown_once():
    """Less than half the GPU free used to be an error; now it is a warning."""
    free, total = mem_info()
    grab = free - total // 2 + 256 * 2**20
    if grab <= 0 or grab > free - 512 * 2**20:
        pytest.skip("cannot leave the GPU between half full and nearly full")
    ptr = ctypes.c_void_p()
    if CUDART.cudaMalloc(ctypes.byref(ptr), grab) != 0:
        pytest.skip("could not allocate GPU memory to make the GPU look busy")
    try:
        geo, angles, img = small_problem()
        proj = tigre.Ax(img, geo, angles)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("default")
            for _ in range(3):
                tigre.Atb(proj, geo, angles)
        busy = [w for w in caught if "heavily used" in str(w.message)]
        assert len(busy) == 1, [str(w.message) for w in caught]
        assert issubclass(busy[0].category, RuntimeWarning)
    finally:
        CUDART.cudaFree(ptr)
