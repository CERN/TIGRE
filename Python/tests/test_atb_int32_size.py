"""Atb allocates its result correctly when geo.nVoxel holds 32-bit integers (e.g. np.array([1024, 1024, 1024])
on Windows with NumPy 1.x, whose default integer is 32-bit).

_Atb.pyx multiplied nVoxel's entries as NumPy scalars, so 1024^3 * 4 bytes wrapped to 0: malloc returned a tiny
block and the backprojection wrote 4 GiB into it (access violation). Needs 4 GiB of GPU and host memory.
"""
import numpy as np
import pytest

import tigre

pytestmark = [pytest.mark.gpu, pytest.mark.extended]  # extended: not run by a plain `pytest`, see contributing.md


def _gpu_free_bytes():
    try:
        import torch
        return torch.cuda.mem_get_info()[0]
    except Exception:
        return None


def test_atb_with_int32_nvoxel_at_1024_cubed():
    free = _gpu_free_bytes()
    if free is not None and free < 5 * 2**30:
        pytest.skip("needs ~5 GiB of free GPU memory")
    n = 1024
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([n, n, n], dtype=np.int32))
    geo.nVoxel = np.array([n, n, n], dtype=np.int32)  # as a user would write it on Windows with NumPy 1.x
    geo.nDetector = np.array([8, 8], dtype=np.int32)
    geo.dDetector = np.array([64.0, 64.0])
    geo.sDetector = geo.nDetector * geo.dDetector
    angles = np.array([0.0, np.pi / 2])
    proj = np.ones((2, 8, 8), dtype=np.float32)
    img = tigre.Atb(proj, geo, angles, "FDK")
    assert img.shape == (n, n, n)
    assert img.dtype == np.float32
    assert np.isfinite(img[::64, ::64, ::64]).all()
    assert img[n // 2, n // 2, n // 2] > 0  # the centre voxel sees both projections
