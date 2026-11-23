"""Conversion between TIGRE cone-beam geometries and per-view (source, detector centre, u, v) vectors.

The CPU tests check round trips in both directions and the refusal of vector geometries TIGRE cannot describe. The
GPU tests check the conversion against the projector itself: line integrals of Gaussian blobs, computed exactly from
the pixel positions the vectors give, must match ``Ax`` on the converted geometry. That check does not depend on
TIGRE's angle conventions, so it catches convention errors the round trips cannot.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from tigre.utilities.vector_geometry import from_TIGRE_to_vector, from_vector_to_TIGRE

N_DETECTOR = np.array([40, 48])  # (rows, columns)


def _base_geo(n_detector=N_DETECTOR, d_detector=(1.3, 1.1)):
    return SimpleNamespace(
        mode="cone",
        nVoxel=np.array([32, 32, 32]),
        sVoxel=np.array([64.0, 64.0, 64.0]),
        dVoxel=np.array([2.0, 2.0, 2.0]),
        nDetector=np.array(n_detector),
        dDetector=np.array(d_detector, dtype=np.float64),
        sDetector=np.array(n_detector) * np.array(d_detector),
    )


def _random_tigre_geo(rng, n):
    geo = _base_geo()
    geo.DSO = rng.uniform(250, 350, n)
    geo.DSD = geo.DSO + rng.uniform(150, 250, n)
    geo.angles = np.stack(
        [rng.uniform(0, 2 * np.pi, n), rng.uniform(-0.4, 0.4, n), rng.uniform(-0.4, 0.4, n)], axis=1
    )
    geo.offDetector = rng.uniform(-8, 8, (n, 2))
    geo.rotDetector = rng.uniform(-0.25, 0.25, (n, 3))
    geo.offOrigin = rng.uniform(-4, 4, (n, 3))
    geo.COR = rng.uniform(-3, 3, n)
    return geo


def _random_rotation(rng, max_angle):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    angle = rng.uniform(-max_angle, max_angle)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _facing_detector(rng, source, distance, du, dv, tilt):
    """Detector centre and (u, v) for a detector on the far side of the image, tilted by up to ``tilt``."""
    s = source / np.linalg.norm(source)
    up = np.array([0.0, 0.0, 1.0]) if abs(s[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e_u = np.cross(up, s)
    e_u /= np.linalg.norm(e_u)
    e_up = np.cross(s, e_u)
    frame = _random_rotation(rng, tilt) @ np.stack([s, e_u, e_up], axis=1)
    centre = -s * distance + rng.uniform(-8, 8) * e_u + rng.uniform(-8, 8) * e_up
    return centre, du * frame[:, 1], dv * frame[:, 2]


def _random_vectors(rng, n, du=1.1, dv=1.3):
    vectors = np.zeros((n, 12))
    for i in range(n):
        direction = rng.normal(size=3)
        direction[2] *= 0.3
        source = direction / np.linalg.norm(direction) * rng.uniform(250, 350)
        centre, u, v = _facing_detector(rng, source, rng.uniform(150, 250), du, dv, 0.3)
        vectors[i] = np.concatenate([source, centre, u, v])
    return vectors


def _static_source_array_vectors(rng, n, du=1.1, dv=1.3):
    """Static linear source array facing one fixed tilted detector: no common rotation centre."""
    source = np.stack(
        [np.full(n, 300.0), np.linspace(-120, 120, n), np.linspace(-10, 10, n)], axis=1
    )
    centre, u, v = _facing_detector(rng, np.array([300.0, 0.0, 0.0]), 200.0, du, dv, 0.2)
    return np.concatenate([source, np.tile(centre, (n, 1)), np.tile(u, (n, 1)), np.tile(v, (n, 1))], axis=1)


# --------------------------------------------------------------------------- round trips


@pytest.mark.parametrize("seed", range(5))
def test_tigre_to_vector_to_tigre_describes_same_geometry(seed):
    rng = np.random.default_rng(seed)
    geo = _random_tigre_geo(rng, 25)
    vectors = from_TIGRE_to_vector(geo)
    geo2 = from_vector_to_TIGRE(geo, vectors)
    np.testing.assert_allclose(from_TIGRE_to_vector(geo2), vectors, rtol=0, atol=1e-9)
    np.testing.assert_allclose(geo2.dDetector, geo.dDetector, rtol=1e-12)


@pytest.mark.parametrize("make", [_random_vectors, _static_source_array_vectors])
@pytest.mark.parametrize("seed", range(3))
def test_vector_to_tigre_to_vector_is_identity(make, seed):
    rng = np.random.default_rng(seed)
    vectors = make(rng, 30)
    geo = from_vector_to_TIGRE(_base_geo(), vectors)
    assert geo.angles.shape == (30, 3)
    assert geo.offOrigin.shape == (30, 3) and not geo.offOrigin.any()
    assert not geo.COR.any()
    np.testing.assert_allclose(geo.dDetector, [1.3, 1.1], rtol=1e-12)
    np.testing.assert_allclose(geo.sDetector, N_DETECTOR * np.array([1.3, 1.1]), rtol=1e-12)
    np.testing.assert_allclose(from_TIGRE_to_vector(geo), vectors, rtol=0, atol=1e-9)


@pytest.mark.parametrize("seed", range(3))
def test_first_pixel_reference(seed):
    rng = np.random.default_rng(seed)
    geo = _random_tigre_geo(rng, 10)
    centre = from_TIGRE_to_vector(geo)
    first = from_TIGRE_to_vector(geo, detector_reference="first_pixel")
    u, v = centre[:, 6:9], centre[:, 9:12]
    expected = centre[:, 3:6] - u * (N_DETECTOR[1] - 1) / 2 - v * (N_DETECTOR[0] - 1) / 2
    np.testing.assert_allclose(first[:, 3:6], expected, atol=1e-9)
    np.testing.assert_allclose(first[:, [0, 1, 2, 6, 7, 8, 9, 10, 11]], centre[:, [0, 1, 2, 6, 7, 8, 9, 10, 11]])
    geo2 = from_vector_to_TIGRE(geo, first, detector_reference="first_pixel")
    np.testing.assert_allclose(from_TIGRE_to_vector(geo2), centre, atol=1e-9)


def test_accepts_column_shaped_distances():
    # Other helpers (e.g. ArbitrarySourceDetectorFixedObject) return DSO/DSD shaped (N, 1).
    geo = _random_tigre_geo(np.random.default_rng(0), 6)
    expected = from_TIGRE_to_vector(geo)
    geo.DSO, geo.DSD = geo.DSO[:, None], geo.DSD[:, None]
    np.testing.assert_array_equal(from_TIGRE_to_vector(geo), expected)


def test_does_not_modify_inputs():
    rng = np.random.default_rng(0)
    geo = _random_tigre_geo(rng, 5)
    before = {k: np.copy(v) for k, v in vars(geo).items() if isinstance(v, np.ndarray)}
    vectors = from_TIGRE_to_vector(geo)
    vectors_before = vectors.copy()
    from_vector_to_TIGRE(geo, vectors)
    for k, v in before.items():
        np.testing.assert_array_equal(getattr(geo, k), v)
    np.testing.assert_array_equal(vectors, vectors_before)


# --------------------------------------------------------------------------- known cases


def test_circular_scan_values():
    geo = _base_geo()
    geo.DSO, geo.DSD = 300.0, 500.0
    angles = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    vectors = from_TIGRE_to_vector(geo, angles)
    c, s = np.cos(angles), np.sin(angles)
    z = np.zeros_like(angles)
    expected = np.stack(
        [300 * c, 300 * s, z, -200 * c, -200 * s, z, -1.1 * s, 1.1 * c, z, z, z, 1.3 + z], axis=1
    )
    np.testing.assert_allclose(vectors, expected, atol=1e-12)

    # A source in the x-y plane comes back as a rotation about z only, with no detector offset or rotation.
    geo2 = from_vector_to_TIGRE(geo, vectors)
    np.testing.assert_allclose(np.unwrap(geo2.angles[:, 0]), angles, atol=1e-12)
    np.testing.assert_allclose(geo2.angles[:, 1:], 0, atol=1e-12)
    np.testing.assert_allclose(geo2.rotDetector, 0, atol=1e-12)
    np.testing.assert_allclose(geo2.offDetector, 0, atol=1e-12)
    np.testing.assert_allclose(geo2.DSO, 300)
    np.testing.assert_allclose(geo2.DSD, 500)


def test_detector_tilt_goes_to_rotDetector():
    geo = _base_geo()
    geo.DSO, geo.DSD = 300.0, 500.0
    geo.offDetector = np.array([3.0, -5.0])
    geo.rotDetector = np.array([0.05, -0.1, 0.2])
    angles = np.linspace(0, np.pi, 7)
    geo2 = from_vector_to_TIGRE(geo, from_TIGRE_to_vector(geo, angles))
    np.testing.assert_allclose(geo2.rotDetector, np.tile(geo.rotDetector, (7, 1)), atol=1e-12)
    np.testing.assert_allclose(geo2.offDetector, np.tile(geo.offDetector, (7, 1)), atol=1e-12)
    np.testing.assert_allclose(geo2.angles[:, 1:], 0, atol=1e-12)


@pytest.mark.parametrize(
    "source",
    [[0.0, 0.0, 300.0], [0.0, 0.0, -300.0]],
    ids=["source above", "source below"],
)
def test_source_on_z_axis(source):
    # The "angles" frame is degenerate (theta = 0 or pi) and its z axis cannot follow world z.
    rng = np.random.default_rng(1)
    source = np.array(source)
    centre, u, v = _facing_detector(rng, source, 200.0, 1.1, 1.3, 0.2)
    vectors = np.concatenate([source, centre, u, v])[None, :]
    geo = from_vector_to_TIGRE(_base_geo(), vectors)
    np.testing.assert_allclose(from_TIGRE_to_vector(geo), vectors, atol=1e-9)


def test_matches_matlab_fixture():
    """Same case and values as MATLAB/Tests/test_vectorGeometry.m, so both languages use one convention."""
    geo = _base_geo()
    geo.DSO, geo.DSD, geo.COR = 300.0, 480.0, 1.5
    geo.offOrigin = np.array([3.0, -2.0, 1.0])  # (z, y, x)
    geo.offDetector = np.array([4.0, -6.0])  # (v, u)
    geo.rotDetector = np.array([0.1, -0.2, 0.3])
    angles = np.array([[0.7, 0.2, -0.3], [2.5, -0.1, 0.4]])
    expected = [
        [2.699824759324622e02, 1.162923904853032e02, -5.993881829360623e01, -1.621642623444666e02,
         -6.973146434746344e01, 3.543582209769166e01, -7.050992103620070e-01, 8.371232783412641e-01,
         1.098167583115110e-01, 8.360878861563337e-02, -9.926385137791681e-02, 1.293505414861434e00],
        [-2.920792266289853e02, 7.174692812602419e01, 2.458579979142952e01, 1.739759740279454e02,
         -3.619112849413428e01, -1.533820143203991e01, 9.139008720326046e-02, -1.095637089932312e00,
         3.503168745057018e-02, 3.499071258227064e-01, 6.915757727574623e-02, 1.250112888024055e00],
    ]
    np.testing.assert_allclose(from_TIGRE_to_vector(geo, angles), expected, rtol=0, atol=1e-9)


# --------------------------------------------------------------------------- refusals


def test_refuses_non_perpendicular_axes():
    vectors = _random_vectors(np.random.default_rng(0), 4)
    u, v = vectors[2, 6:9], vectors[2, 9:12]
    vectors[2, 9:12] = np.cos(0.01) * v + np.sin(0.01) * np.linalg.norm(v) / np.linalg.norm(u) * u  # same length
    with pytest.raises(ValueError, match="not perpendicular in view 2"):
        from_vector_to_TIGRE(_base_geo(), vectors)


@pytest.mark.parametrize("column", [6, 9])
def test_refuses_changing_pixel_size(column):
    vectors = _random_vectors(np.random.default_rng(0), 4)
    vectors[3, column : column + 3] *= 1.001
    with pytest.raises(ValueError, match="pixel size"):
        from_vector_to_TIGRE(_base_geo(), vectors)


def test_refuses_source_at_image_centre():
    vectors = _random_vectors(np.random.default_rng(0), 3)
    vectors[1, 0:3] = 0
    with pytest.raises(ValueError, match="centre of the image"):
        from_vector_to_TIGRE(_base_geo(), vectors)


def test_refuses_source_in_detector_plane():
    vectors = _random_vectors(np.random.default_rng(0), 3)
    vectors[0, 0:3] = vectors[0, 3:6] + 50 * vectors[0, 6:9]
    with pytest.raises(ValueError, match="plane of the detector"):
        from_vector_to_TIGRE(_base_geo(), vectors)


@pytest.mark.parametrize(
    "vectors",
    [np.zeros((3, 11)), np.zeros((0, 12)), np.full((2, 12), np.nan)],
    ids=["wrong width", "empty", "nan"],
)
def test_refuses_bad_arrays(vectors):
    with pytest.raises(ValueError):
        from_vector_to_TIGRE(_base_geo(), vectors)


def test_refuses_bad_reference_and_parallel_mode():
    geo = _random_tigre_geo(np.random.default_rng(0), 2)
    with pytest.raises(ValueError, match="detector_reference"):
        from_TIGRE_to_vector(geo, detector_reference="corner")
    geo.mode = "parallel"
    with pytest.raises(ValueError, match="cone-beam"):
        from_TIGRE_to_vector(geo)


# --------------------------------------------------------------------------- projector check (GPU)

# Gaussian blobs (centre x, y, z in mm, sigma in mm, peak value), placed off-centre so a mirrored or swapped axis
# changes the result. Smooth, so voxelising them costs little accuracy, and their line integrals are exact.
BLOBS = [((10.0, -6.0, 4.0), 7.0, 1.0), ((-14.0, 12.0, -9.0), 5.0, 2.0), ((4.0, 16.0, 15.0), 4.0, 3.0)]


def _phantom(geo):
    """BLOBS sampled at the voxel centres. Image array order is (z, y, x)."""
    n = np.asarray(geo.nVoxel)
    d = np.asarray(geo.dVoxel, dtype=np.float64)
    s = np.asarray(geo.sVoxel, dtype=np.float64)
    axes = [-s[k] / 2 + d[k] * (np.arange(n[k]) + 0.5) for k in range(3)]
    z, y, x = np.meshgrid(*axes, indexing="ij")
    img = np.zeros(z.shape)
    for (cx, cy, cz), sigma, value in BLOBS:
        img += value * np.exp(-((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) / (2 * sigma**2))
    return img.astype(np.float32)


def _exact_projections(vectors, n_detector):
    """Line integrals of BLOBS from each source to each pixel centre the vectors define."""
    rows, columns = np.arange(n_detector[0]), np.arange(n_detector[1])
    out = np.zeros((len(vectors), n_detector[0], n_detector[1]))
    for i, vec in enumerate(vectors):
        source, centre, u, v = vec[0:3], vec[3:6], vec[6:9], vec[9:12]
        first = centre - u * (n_detector[1] - 1) / 2 - v * (n_detector[0] - 1) / 2
        pixels = first + rows[:, None, None] * v + columns[None, :, None] * u
        ray = pixels - source
        ray /= np.linalg.norm(ray, axis=-1, keepdims=True)
        for c, sigma, value in BLOBS:
            w = np.asarray(c) - source
            distance2 = w @ w - np.sum(ray * w, axis=-1) ** 2  # squared distance from the blob centre to the ray
            out[i] += value * sigma * np.sqrt(2 * np.pi) * np.exp(-distance2 / (2 * sigma**2))
    return out


def _tigre_geometry_for_gpu(geo):
    import tigre

    full = tigre.geometry(mode="cone", default=True, nVoxel=np.array([64, 64, 64]))
    full.sVoxel = np.array([80.0, 80.0, 80.0])
    full.dVoxel = full.sVoxel / full.nVoxel
    full.nDetector = np.array([60, 72])
    full.dDetector = np.array([1.3, 1.1])
    full.sDetector = full.nDetector * full.dDetector
    full.accuracy = 0.5
    for name in ("DSO", "DSD", "offOrigin", "offDetector", "rotDetector", "COR"):
        if hasattr(geo, name):
            setattr(full, name, getattr(geo, name))
    return full


def _check_projections(geo, angles, vectors):
    import tigre

    proj = tigre.Ax(_phantom(geo), geo, angles, projection_type="interpolated")
    exact = _exact_projections(vectors, geo.nDetector)
    # Sampling the blobs on the voxel grid and the projector discretisation leave about 0.8%; a detector off by half
    # a pixel gives 4-5%, and a wrong sign or axis tens of percent or more.
    error = np.linalg.norm(proj - exact) / np.linalg.norm(exact)
    assert error < 0.015, error
    return proj


@pytest.mark.gpu
def test_vectors_from_tigre_match_projector():
    rng = np.random.default_rng(7)
    geo = _tigre_geometry_for_gpu(_random_tigre_geo(rng, 12))
    angles = _random_tigre_geo(rng, 12).angles
    # Image is offset by offOrigin, so keep the spheres in image-centred coordinates (as the vectors are).
    _check_projections(geo, angles, from_TIGRE_to_vector(geo, angles))


@pytest.mark.gpu
@pytest.mark.parametrize("make", [_random_vectors, _static_source_array_vectors])
def test_tigre_from_vectors_matches_projector(make):
    vectors = make(np.random.default_rng(3), 12)
    geo = from_vector_to_TIGRE(_tigre_geometry_for_gpu(SimpleNamespace()), vectors)
    _check_projections(geo, geo.angles, vectors)


@pytest.mark.gpu
def test_equivalent_geometries_give_same_operators():
    """The original geometry and its round trip are different TIGRE descriptions of the same scan.

    COR is zero here: the backprojectors apply COR along the rotated detector's u axis, while the projectors apply
    it along (-sin(alpha), cos(alpha), 0), so with a rotated detector or nonzero theta/psi they disagree whatever the
    conversion does. The "FDK" backprojector is not compared: its weights depend on DSO and the rotation axis, so
    equivalent descriptions of one scan legitimately weight it differently.
    """
    import tigre

    rng = np.random.default_rng(11)
    geo = _tigre_geometry_for_gpu(_random_tigre_geo(rng, 10))
    geo.COR = 0.0
    angles = _random_tigre_geo(rng, 10).angles
    geo2 = from_vector_to_TIGRE(geo, from_TIGRE_to_vector(geo, angles))
    img = _phantom(geo)
    p1, p2 = tigre.Ax(img, geo, angles), tigre.Ax(img, geo2, geo2.angles)
    assert np.linalg.norm(p1 - p2) / np.linalg.norm(p1) < 1e-4
    b1 = tigre.Atb(p1, geo, angles, backprojection_type="matched")
    b2 = tigre.Atb(p1, geo2, geo2.angles, backprojection_type="matched")
    assert np.linalg.norm(b1 - b2) / np.linalg.norm(b1) < 1e-4
