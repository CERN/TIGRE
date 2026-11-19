"""Atb must not leak texture objects (https://github.com/CERN/TIGRE/issues/617).

Atb uploads the projections in blocks of 32 angles and creates a texture object
for each block in one of two slots, but only destroyed the two left at the end.
Every call leaked (blocks - 2) texture objects, until cudaCreateTextureObject
failed with "out of memory" after about a million of them: around iteration 707
of the reproduction in the issue, regardless of how much memory was free.

Each texture object also holds a small amount of device memory (~60 bytes in
the reproduction), so the leak shows up as a steady loss of free memory. With
many angles and many calls the old code loses several MB; the fixed code none.
"""
import ctypes
import os

import numpy as np
import pytest

import tigre

pytestmark = [pytest.mark.gpu, pytest.mark.extended]  # extended: not run by a plain `pytest`, see contributing.md

BLOCK = 32                 # PROJ_PER_KERNEL in the backprojectors
NANGLES = BLOCK * 400      # 400 blocks: 398 texture objects leaked per call before the fix
NCALLS = 500               # ~200k leaked objects, roughly 6-12 MB, before the fix
MAX_DRIFT_MB = 3.0


def _cudart():
    names = ["cudart64_13.dll", "cudart64_12.dll", "cudart64_110.dll",
             "libcudart.so", "libcudart.so.13", "libcudart.so.12", "libcudart.so.11.0"]
    dirs = [""]
    if os.environ.get("CUDA_PATH"):
        dirs.append(os.path.join(os.environ["CUDA_PATH"], "bin"))
    for d in dirs:
        for n in names:
            try:
                return ctypes.CDLL(os.path.join(d, n) if d else n)
            except OSError:
                continue
    return None


_CUDART = _cudart()


def gpu_free_bytes():
    try:
        import torch
        if torch.cuda.is_available():
            return torch.cuda.mem_get_info()[0]
    except Exception:  # noqa: BLE001  torch missing or broken: use the runtime directly
        pass
    if _CUDART is None:
        return None
    free, total = ctypes.c_size_t(0), ctypes.c_size_t(0)
    if _CUDART.cudaMemGetInfo(ctypes.byref(free), ctypes.byref(total)) != 0:
        return None
    return free.value


@pytest.mark.parametrize("mode,backprojection_type", [
    ("cone", "matched"),
    ("cone", "FDK"),
    ("parallel", "FDK"),
])
def test_atb_does_not_leak_texture_objects(mode, backprojection_type):
    if gpu_free_bytes() is None:
        pytest.skip("cannot read free GPU memory (no torch and CUDA runtime not loadable)")

    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array([8, 8, 8]))
    angles = np.linspace(0, 2 * np.pi, NANGLES, endpoint=False, dtype=np.float32)
    proj = np.ones((NANGLES,) + tuple(geo.nDetector), dtype=np.float32)

    tigre.Atb(proj, geo, angles, backprojection_type)  # absorb one-off setup
    before = gpu_free_bytes()
    for _ in range(NCALLS):
        tigre.Atb(proj, geo, angles, backprojection_type)
    after = gpu_free_bytes()

    drift_mb = (before - after) / 2**20
    assert drift_mb < MAX_DRIFT_MB, (
        "free GPU memory dropped {:.1f} MB over {} Atb calls with {} angles: "
        "texture objects are leaking".format(drift_mb, NCALLS, NANGLES))
