"""COR means the same thing in every projector and backprojector.

``geo.COR`` moves the whole scan (source and detector) by ``COR * (-sin(alpha), cos(alpha), 0)``, so it is the same as
moving the image the other way with ``offOrigin``. The projectors did this; the backprojectors applied COR along the
rotated detector's u axis instead, which only agrees for an unrotated detector and theta = psi = 0, and the FDK weight
added COR to the voxel's y coordinate. Each operator must now give the same result for both descriptions.
"""
import numpy as np
import pytest

import tigre

pytestmark = pytest.mark.gpu

N_ANGLES = 24
TOL = 1e-4  # float32 sums in a different order; any inconsistency is percent-level


def _geometry(mode, rng):
    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array([48, 48, 48]))
    geo.sVoxel = np.array([96.0, 96.0, 96.0])
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    geo.nDetector = np.array([56, 72])
    geo.dDetector = np.array([1.6, 1.6]) if mode == "cone" else np.array([1.2, 1.2])
    geo.sDetector = geo.nDetector * geo.dDetector
    # The interpolated projector samples each ray in steps from where it starts, which differs between the two
    # descriptions; finer steps keep that difference (2e-4 at the default 0.5, parallel beam) well below TOL.
    geo.accuracy = 0.05
    geo.offOrigin = rng.uniform(-3, 3, (N_ANGLES, 3))
    geo.COR = rng.uniform(-4, 4, N_ANGLES)
    geo.rotDetector = rng.uniform(-0.15, 0.15, (N_ANGLES, 3))
    if mode == "cone":
        geo.DSO, geo.DSD = 300.0, 480.0
        geo.offDetector = rng.uniform(-5, 5, (N_ANGLES, 2))
    return geo


def _without_cor(geo, angles):
    """The same scan with COR = 0 and the image moved by -COR * (-sin(alpha), cos(alpha), 0) instead."""
    moved = tigre.geometry(mode=geo.mode, default=True, nVoxel=geo.nVoxel)
    for name, value in vars(geo).items():
        setattr(moved, name, np.copy(value) if isinstance(value, np.ndarray) else value)
    alpha = angles[:, 0]
    offset = np.array(geo.offOrigin, dtype=np.float64)  # (z, y, x)
    offset[:, 2] += geo.COR * np.sin(alpha)
    offset[:, 1] -= geo.COR * np.cos(alpha)
    moved.offOrigin = offset
    moved.COR = np.zeros(N_ANGLES)
    return moved


def _angles(mode, rng):
    alpha = np.linspace(0, 2 * np.pi, N_ANGLES, endpoint=False)
    return np.stack([alpha, rng.uniform(-0.2, 0.2, N_ANGLES), rng.uniform(-0.2, 0.2, N_ANGLES)], axis=1)


def _image(geo):
    z, y, x = np.meshgrid(*[np.linspace(-1, 1, n) for n in geo.nVoxel], indexing="ij")
    img = np.exp(-((x - 0.2) ** 2 + (y + 0.3) ** 2 + (z - 0.1) ** 2) / 0.1)
    img += 2 * np.exp(-((x + 0.4) ** 2 + (y - 0.2) ** 2 + (z + 0.3) ** 2) / 0.03)
    return img.astype(np.float32)


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(a)


@pytest.mark.parametrize("mode", ["cone", "parallel"])
@pytest.mark.parametrize("projection_type", ["interpolated", "Siddon"])
def test_projectors(mode, projection_type):
    rng = np.random.default_rng(0)
    geo, angles = _geometry(mode, rng), _angles(mode, rng)
    img = _image(geo)
    a = tigre.Ax(img, geo, angles, projection_type=projection_type)
    b = tigre.Ax(img, _without_cor(geo, angles), angles, projection_type=projection_type)
    assert _rel(a, b) < TOL


@pytest.mark.parametrize(
    "mode, backprojection_type", [("cone", "FDK"), ("cone", "matched"), ("parallel", "FDK"), ("parallel", "matched")]
)
def test_backprojectors(mode, backprojection_type):
    rng = np.random.default_rng(1)
    geo, angles = _geometry(mode, rng), _angles(mode, rng)
    proj = tigre.Ax(_image(geo), geo, angles)
    a = tigre.Atb(proj, geo, angles, backprojection_type=backprojection_type)
    b = tigre.Atb(proj, _without_cor(geo, angles), angles, backprojection_type=backprojection_type)
    assert _rel(a, b) < TOL
