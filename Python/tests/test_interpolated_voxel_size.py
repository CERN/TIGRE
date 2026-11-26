"""The interpolated projector agrees with Siddon for non-integer voxel sizes.

With theta or psi nonzero, the interpolated projector skips the start of each
ray by a distance computed on the host from the largest voxel size. If that
size is rounded to an integer, sampling starts too far along the ray: voxel
sizes above 1 drop part or all of the object, and sizes below 1 divide by
zero.

Needs a GPU. 64^3 at 12 views.
"""
import numpy as np
import pytest

import tigre

N = 64


@pytest.fixture(scope="module")
def phantom():
    z, y, x = np.mgrid[:N, :N, :N].astype(np.float32)
    img = np.zeros((N, N, N), dtype=np.float32)
    img[((x - 30) ** 2 + (y - 32) ** 2 + (z - 32) ** 2) < 14 ** 2] = 1.0
    img[((x - 40) ** 2 + (y - 38) ** 2 + (z - 30) ** 2) < 6 ** 2] = 2.0
    return img


@pytest.mark.parametrize("mode", ["cone", "parallel"])
@pytest.mark.parametrize("dvoxel", [0.8, 1.5])
def test_interpolated_matches_siddon_with_euler_angles(phantom, mode, dvoxel):
    geo = tigre.geometry(mode=mode, nVoxel=np.array([N, N, N]), default=True)
    geo.DSO = np.float32(1000)
    geo.DSD = np.float32(1536)
    geo.dVoxel = np.array([dvoxel, dvoxel, dvoxel])
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.nDetector = np.array([2 * N, 2 * N])
    geo.dDetector = np.array([1.6, 1.6]) * dvoxel
    geo.sDetector = geo.nDetector * geo.dDetector
    n = 12
    angles = np.stack([np.linspace(0, 2 * np.pi, n, endpoint=False),
                       np.full(n, 0.3), np.full(n, 0.2)], axis=1).astype(np.float32)

    siddon = tigre.Ax(phantom, geo, angles, "Siddon")
    interp = tigre.Ax(phantom, geo, angles, "interpolated")

    assert np.isfinite(interp).all()
    err = float(np.linalg.norm(interp - siddon) / np.linalg.norm(siddon))
    assert err < 0.05, f"interpolated vs Siddon relative error {err:.3f}"
