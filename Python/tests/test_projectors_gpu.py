"""GPU tests for the forward projectors (tigre.Ax) and backprojectors (tigre.Atb).

TIGRE's projector pairs are not exact adjoints of each other (the default "FDK"
backprojector is not meant to be, and "matched" is only approximately matched),
so these tests do not check <Ax, y> == <x, A^T y>. They check properties every
pair must have: shapes, linearity, non-negativity, determinism, known line
integrals, agreement between the two forward projectors, and that splitting the
work over several GPUs does not change the result.

Run them locally on a machine with a GPU:  pytest -m gpu Python/tests
"""
import numpy as np
import pytest

import tigre
from tigre.utilities import sample_loader
from tigre.utilities.gpu import GpuIds

pytestmark = pytest.mark.gpu

N = 64
N_ANGLES = 36

FORWARD_TYPES = ["Siddon", "interpolated"]
BACKWARD_TYPES = ["FDK", "matched"]


def cone_geometry():
    # Low-resolution default: 64^3 volume of 256 mm, 128x128 detector that sees all of it.
    return tigre.geometry(mode="cone", default=True, high_resolution=False)


def parallel_geometry():
    # 64^3 volume with 1 mm voxels and a 64x64 detector with 1 mm pixels.
    return tigre.geometry(mode="parallel", nVoxel=np.array([N, N, N]))


GEOMETRIES = {"cone": cone_geometry, "parallel": parallel_geometry}


def angles(n=N_ANGLES):
    return np.linspace(0, 2 * np.pi, n, endpoint=False, dtype=np.float32)


def head(geo):
    return sample_loader.load_head_phantom(geo.nVoxel).astype(np.float32)


def random_volume(geo, seed):
    return np.random.default_rng(seed).random(tuple(geo.nVoxel), dtype=np.float32)


def random_projections(geo, n_angles, seed):
    shape = (n_angles, int(geo.nDetector[0]), int(geo.nDetector[1]))
    return np.random.default_rng(seed).random(shape, dtype=np.float32)


def rel_diff(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


@pytest.fixture(params=sorted(GEOMETRIES))
def geo(request):
    return GEOMETRIES[request.param]()


# --- shapes and trivial inputs -----------------------------------------------------------


@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_shape_and_dtype(geo, ptype):
    a = angles()
    proj = tigre.Ax(head(geo), geo, a, ptype)
    assert proj.shape == (len(a), geo.nDetector[0], geo.nDetector[1])
    assert proj.dtype == np.float32
    assert np.isfinite(proj).all()


@pytest.mark.parametrize("btype", BACKWARD_TYPES)
def test_backward_shape_and_dtype(geo, btype):
    a = angles()
    vol = tigre.Atb(random_projections(geo, len(a), 0), geo, a, btype)
    assert vol.shape == tuple(geo.nVoxel)
    assert vol.dtype == np.float32
    assert np.isfinite(vol).all()


@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_of_zero_is_zero(geo, ptype):
    a = angles()
    proj = tigre.Ax(np.zeros(tuple(geo.nVoxel), np.float32), geo, a, ptype)
    assert not proj.any()


@pytest.mark.parametrize("btype", BACKWARD_TYPES)
def test_backward_of_zero_is_zero(geo, btype):
    a = angles()
    shape = (len(a), geo.nDetector[0], geo.nDetector[1])
    vol = tigre.Atb(np.zeros(shape, np.float32), geo, a, btype)
    assert not vol.any()


# --- linearity and positivity ------------------------------------------------------------


@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_is_linear(geo, ptype):
    a = angles()
    x, y = random_volume(geo, 1), random_volume(geo, 2)
    lhs = tigre.Ax(2 * x + 3 * y, geo, a, ptype)
    rhs = 2 * tigre.Ax(x, geo, a, ptype) + 3 * tigre.Ax(y, geo, a, ptype)
    assert rel_diff(lhs, rhs) < 1e-5


@pytest.mark.parametrize("btype", BACKWARD_TYPES)
def test_backward_is_linear(geo, btype):
    a = angles()
    p, q = random_projections(geo, len(a), 3), random_projections(geo, len(a), 4)
    lhs = tigre.Atb(2 * p + 3 * q, geo, a, btype)
    rhs = 2 * tigre.Atb(p, geo, a, btype) + 3 * tigre.Atb(q, geo, a, btype)
    assert rel_diff(lhs, rhs) < 1e-5


@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_preserves_nonnegativity(geo, ptype):
    proj = tigre.Ax(random_volume(geo, 5), geo, angles(), ptype)
    assert proj.min() >= 0
    assert proj.max() > 0


@pytest.mark.parametrize("btype", BACKWARD_TYPES)
def test_backward_preserves_nonnegativity(geo, btype):
    a = angles()
    vol = tigre.Atb(random_projections(geo, len(a), 6), geo, a, btype)
    assert vol.min() >= 0
    assert vol.max() > 0


# --- determinism -------------------------------------------------------------------------


@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_is_deterministic(geo, ptype):
    a = angles()
    x = head(geo)
    np.testing.assert_array_equal(tigre.Ax(x, geo, a, ptype), tigre.Ax(x, geo, a, ptype))


@pytest.mark.parametrize("btype", BACKWARD_TYPES)
def test_backward_is_deterministic(geo, btype):
    a = angles()
    p = random_projections(geo, len(a), 7)
    np.testing.assert_array_equal(tigre.Atb(p, geo, a, btype), tigre.Atb(p, geo, a, btype))


# --- known values ------------------------------------------------------------------------


@pytest.mark.parametrize("ptype,rtol", [("Siddon", 1e-4), ("interpolated", 2e-2)])
@pytest.mark.parametrize("angle", [0.0, np.pi / 2, np.pi, 3 * np.pi / 2])
def test_parallel_line_integral_through_uniform_cube(ptype, rtol, angle):
    """With 1 mm voxels, every central ray crosses N voxels of value 1 along an axis."""
    geo = parallel_geometry()
    a = np.array([angle], dtype=np.float32)
    proj = tigre.Ax(np.ones((N, N, N), np.float32), geo, a, ptype)
    centre = proj[0, N // 2 - 4 : N // 2 + 4, N // 2 - 4 : N // 2 + 4]
    np.testing.assert_allclose(centre, N, rtol=rtol)


def test_siddon_and_interpolated_agree(geo):
    a = angles()
    x = head(geo)
    siddon = tigre.Ax(x, geo, a, "Siddon")
    interp = tigre.Ax(x, geo, a, "interpolated")
    assert rel_diff(interp, siddon) < 0.05


@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_is_invariant_to_a_full_turn(geo, ptype):
    x = head(geo)
    a = np.array([0.3, 1.1], dtype=np.float32)
    first = tigre.Ax(x, geo, a, ptype)
    turned = tigre.Ax(x, geo, a + np.float32(2 * np.pi), ptype)
    assert rel_diff(turned, first) < 1e-4


# --- multi-GPU ---------------------------------------------------------------------------


def _two_or_more_gpus():
    try:
        return len(GpuIds()) >= 2
    except Exception:
        return False


multi_gpu = pytest.mark.skipif(not _two_or_more_gpus(), reason="needs at least 2 GPUs")


@multi_gpu
@pytest.mark.parametrize("ptype", FORWARD_TYPES)
def test_forward_multi_gpu_matches_single_gpu(geo, ptype):
    a = angles()
    x = head(geo)
    single = tigre.Ax(x, geo, a, ptype, gpuids=GpuIds()[0])
    multi = tigre.Ax(x, geo, a, ptype, gpuids=GpuIds())
    np.testing.assert_allclose(multi, single, rtol=1e-6, atol=1e-6 * np.abs(single).max())


@multi_gpu
@pytest.mark.parametrize("btype", BACKWARD_TYPES)
def test_backward_multi_gpu_matches_single_gpu(geo, btype):
    a = angles()
    p = random_projections(geo, len(a), 8)
    single = tigre.Atb(p, geo, a, btype, gpuids=GpuIds()[0])
    multi = tigre.Atb(p, geo, a, btype, gpuids=GpuIds())
    np.testing.assert_allclose(multi, single, rtol=1e-6, atol=1e-6 * np.abs(single).max())
