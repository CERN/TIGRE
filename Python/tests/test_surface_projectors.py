"""Native surface rays checked against independent float64 voxel intersections."""

from types import SimpleNamespace

import numpy as np
import pytest

native = pytest.importorskip("_AxSurface")
import tigre
from tigre.utilities.Ax_surface import Ax_surface, Atb_surface
from tigre.utilities.detector_geometry import DetectorSurfaceGeometry
from tigre.utilities.gpu import getGpuNames

pytestmark = pytest.mark.skipif(not getGpuNames(), reason="requires a CUDA GPU")


def fixture(curved=True):
    geo = SimpleNamespace(
        nVoxel=np.array([4, 5, 6]),
        dVoxel=np.array([1.5, 0.8, 1.2]),
        offOrigin=np.array([0.17, -0.21, 0.13]),
    )
    beta = np.array([0.1, 0.7, 1.6])
    rot = np.zeros((3, 3, 3))
    rot[:, 0, 0] = rot[:, 1, 1] = np.cos(beta)
    rot[:, 1, 0] = np.sin(beta)
    rot[:, 0, 1] = -np.sin(beta)
    rot[:, 2, 2] = 1.0
    source = rot @ np.array([30.0, 0.0, 0.3])
    local = np.zeros((3, 5, 3))
    u = np.linspace(-7, 7, 5)
    local[..., 0] = -20.0 + (0.03 * u * u if curved else 0)
    local[..., 1] = u
    local[..., 2] = np.linspace(-4.0, 4.0, 3)[:, None]
    if curved == "u":
        local[..., 0] = [-7.0, -19.0, -20.0, -19.0, -7.0]
        local[..., 1] = [-4.0, -4.0, 0.0, 4.0, 4.0]
    pixels = np.einsum("vij,rcj->vrci", rot, local)
    mask = np.ones((3, 3, 5), bool)
    mask[:, :, 2] = False
    return geo, DetectorSurfaceGeometry(source, pixels, mask)


def dense_ray_matrix(geo, surface):
    """Intersect every voxel box separately; no production traversal is reused."""
    shape = tuple(geo.nVoxel)
    spacing = geo.dVoxel[::-1]
    origin = geo.offOrigin[::-1] - 0.5 * geo.nVoxel[::-1] * spacing
    matrix = np.zeros((surface.valid_mask.size, int(np.prod(shape))))
    for flat, index in enumerate(np.ndindex(surface.valid_mask.shape)):
        if not surface.valid_mask[index]:
            continue
        source = surface.source_positions_mm[index[0]].astype(float)
        ray = surface.pixel_centers_mm[index].astype(float) - source
        length = np.linalg.norm(ray)
        for col, voxel in enumerate(np.ndindex(shape)):
            lo = origin + np.array(voxel[::-1]) * spacing
            hi = lo + spacing
            enter, leave = 0.0, 1.0
            for axis in range(3):
                if ray[axis] == 0:
                    if not lo[axis] <= source[axis] < hi[axis]:
                        leave = -1
                        break
                else:
                    a, b = sorted(
                        (
                            (lo[axis] - source[axis]) / ray[axis],
                            (hi[axis] - source[axis]) / ray[axis],
                        )
                    )
                    enter, leave = max(enter, a), min(leave, b)
            matrix[flat, col] = max(0.0, leave - enter) * length
    return matrix


@pytest.mark.parametrize("curved", [False, True, "u"])
def test_forward_and_adjoint_match_independent_dense_matrix(curved):
    geo, surface = fixture(curved)
    matrix = dense_ray_matrix(geo, surface)
    rng = np.random.default_rng(7)
    image = rng.normal(size=geo.nVoxel).astype(np.float32)
    data = rng.normal(size=surface.valid_mask.shape).astype(np.float32)
    actual = Ax_surface(image, geo, surface, gpuids=0)
    back = Atb_surface(data, geo, surface, gpuids=0)
    np.testing.assert_allclose(actual.ravel(), matrix @ image.ravel(), atol=2e-5, rtol=2e-4)
    np.testing.assert_allclose(back.ravel(), matrix.T @ data.ravel(), atol=2e-5, rtol=2e-4)
    left, right = np.vdot(actual.astype(float), data), np.vdot(image.astype(float), back)
    assert abs(left - right) / max(abs(left), abs(right), 1.0) < 2e-5


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_axis_aligned_uniform_box_integral(axis):
    geo, _ = fixture()
    source = geo.offOrigin[::-1].copy()
    pixel = source.copy()
    source[axis], pixel[axis] = 30.0, -30.0
    surface = DetectorSurfaceGeometry(source[None], pixel.reshape(1, 1, 3))
    value = Ax_surface(np.ones(geo.nVoxel, np.float32), geo, surface, gpuids=0)
    assert value.item() == pytest.approx((geo.nVoxel * geo.dVoxel)[::-1][axis], abs=2e-5)


def test_masked_nonfinite_data_and_slab_decomposition():
    geo, surface = fixture()
    image = np.ones(geo.nVoxel, np.float32)
    single = Ax_surface(image, geo, surface, gpuids=0)
    split = Ax_surface(image, geo, surface, gpuids=0, n_slabs=2)
    np.testing.assert_allclose(split, single, atol=2e-5)
    dirty = single.copy()
    dirty[~surface.valid_mask] = np.nan
    before = dirty.copy()
    np.testing.assert_allclose(
        Atb_surface(dirty, geo, surface, gpuids=0),
        Atb_surface(single, geo, surface, gpuids=0),
        atol=2e-5,
    )
    np.testing.assert_array_equal(dirty, before)
    np.testing.assert_allclose(
        Atb_surface(single, geo, surface, gpuids=0, n_slabs=2),
        Atb_surface(single, geo, surface, gpuids=0),
        atol=2e-5,
    )


def test_invalid_device_recovers_and_malformed_shapes_raise():
    geo, surface = fixture()
    image = np.ones(geo.nVoxel, np.float32)
    for _ in range(3):
        with pytest.raises(ValueError):
            Ax_surface(image, geo, surface, gpuids=99999)
        assert np.isfinite(Ax_surface(image, geo, surface, gpuids=0)).all()
    with pytest.raises(ValueError):
        native._Atb_surface_ext(
            np.zeros((3, 1, 1), np.float32),
            tuple(geo.nVoxel),
            tuple(geo.dVoxel[::-1]),
            tuple(geo.offOrigin[::-1]),
            surface.source_positions_mm,
            surface.pixel_centers_mm,
            None,
            0,
        )


def test_flat_detector_matches_existing_siddon_forward():
    geo = tigre.geometry(default=True, nVoxel=np.array([16, 16, 16]))
    geo.DSO, geo.DSD = 100.0, 200.0
    geo.nVoxel = np.array([16, 16, 16])
    geo.dVoxel = np.ones(3)
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.offOrigin = np.zeros(3)
    geo.nDetector = np.array([20, 24])
    geo.dDetector = np.array([1.5, 1.2])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offDetector = np.zeros(2)
    geo.rotDetector = np.zeros(3)
    geo.COR = 0.0
    beta = np.linspace(0.1, 2 * np.pi + 0.1, 9, endpoint=False).astype(np.float32)
    local = np.zeros((20, 24, 3))
    local[..., 0] = -100.0
    local[..., 1] = (np.arange(24) - 11.5)[None] * 1.2
    local[..., 2] = (np.arange(20) - 9.5)[:, None] * 1.5
    rotation = np.array(
        [[[np.cos(b), -np.sin(b), 0], [np.sin(b), np.cos(b), 0], [0, 0, 1]] for b in beta]
    )
    surface = DetectorSurfaceGeometry(
        rotation @ np.array([100.0, 0.0, 0.0]), np.einsum("vij,rcj->vrci", rotation, local)
    )
    z, y, x = np.meshgrid(*[np.arange(16) - 7.5] * 3, indexing="ij")
    image = np.exp(-((x - 2) ** 2 + (y + 1) ** 2 + (z - 1) ** 2) / 8).astype(np.float32)
    expected = tigre.Ax(image, geo, beta)
    actual = Ax_surface(image, geo, surface, gpuids=0)
    assert np.linalg.norm(actual - expected) / np.linalg.norm(expected) < 2e-3


@pytest.mark.parametrize("devices", [True, 0.5, [-1], [], [2**40]])
def test_invalid_device_selection_is_rejected(devices):
    geo, surface = fixture()
    with pytest.raises(ValueError, match="gpuids"):
        Ax_surface(np.ones(geo.nVoxel, np.float32), geo, surface, gpuids=devices)


@pytest.mark.parametrize("slabs", [True, 0, -1, 1.5])
def test_invalid_slab_count_is_rejected(slabs):
    geo, surface = fixture()
    with pytest.raises(ValueError, match="n_slabs"):
        Ax_surface(np.ones(geo.nVoxel, np.float32), geo, surface, gpuids=0, n_slabs=slabs)


@pytest.mark.parametrize("spacing", [(1e100, 1.0, 1.0), (1e-100, 1.0, 1.0)])
def test_native_spacing_must_remain_positive_finite_float32(spacing):
    geo, surface = fixture()
    with pytest.raises(ValueError, match="voxel_size"):
        native._Ax_surface_ext(
            np.ones(geo.nVoxel, np.float32),
            spacing,
            (0.0, 0.0, 0.0),
            surface.source_positions_mm,
            surface.pixel_centers_mm,
            None,
            0,
        )


@pytest.mark.parametrize(
    "source,pixel",
    [
        ([5.0, 0.0, 0.0], [0.0, 0.0, 0.0]),  # detector inside the volume
        ([0.0, 0.0, 0.0], [-5.0, 0.0, 0.0]),  # source inside the volume
        ([0.2, 0.3, 0.4], [0.4, 0.5, 0.6]),  # entire segment in one voxel
        ([5.0, 5.0, 5.0], [-5.0, -5.0, -5.0]),  # simultaneous voxel-plane crossings
        ([5.0, -2.0, 0.0], [-5.0, -2.0, 0.0]),  # ray on the included lower face
        ([5.0, 2.0, 0.0], [-5.0, 2.0, 0.0]),  # ray on the excluded upper face
        ([5.0, 3.0, 0.0], [-5.0, 3.0, 0.0]),  # parallel ray outside volume
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),  # zero-length segment
    ],
)
def test_segment_clipping_and_boundary_rays_match_oracle(source, pixel):
    geo = SimpleNamespace(nVoxel=np.array([4, 4, 4]), dVoxel=np.ones(3), offOrigin=np.zeros(3))
    surface = DetectorSurfaceGeometry([source], np.asarray(pixel).reshape(1, 1, 3))
    matrix = dense_ray_matrix(geo, surface)
    image = np.arange(64, dtype=np.float32).reshape(4, 4, 4) / 64
    np.testing.assert_allclose(
        Ax_surface(image, geo, surface, gpuids=0).ravel(),
        matrix @ image.ravel(),
        atol=2e-5,
        rtol=2e-4,
    )
    data = np.ones((1, 1, 1), np.float32)
    np.testing.assert_allclose(
        Atb_surface(data, geo, surface, gpuids=0).ravel(),
        matrix.T @ data.ravel(),
        atol=2e-5,
        rtol=2e-4,
    )
