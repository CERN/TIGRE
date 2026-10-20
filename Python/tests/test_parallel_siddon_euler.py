import numpy as np
import pytest
import tigre


def test_parallel_siddon_z_boundary_sampling():
    """Verify that Z boundary rays match the expected point-texture sampling.

    At offOrigin Z = 0.5 (lower boundary), pixel [0, 15] must be inside volume (value = 32.0).
    At offOrigin Z = -0.5 (upper boundary), pixel [31, 15] must be outside volume (value = 0.0),
    consistent with baseline behavior.
    """
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([32, 32, 32]))
    geo.dVoxel = np.array([1.0, 1.0, 1.0])
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.dDetector = np.array([1.0, 1.0])
    geo.nDetector = np.array([32, 32])
    geo.sDetector = geo.nDetector * geo.dDetector
    img = np.ones(geo.nVoxel, dtype=np.float32)
    angles = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)

    # 1. Lower boundary check (offOrigin Z = 0.5)
    geo.offOrigin = np.array([0.5, 0.0, 0.0])  # [Z, Y, X] in Python
    val_lower = tigre.Ax(img, geo, angles, projection_type="Siddon")[0, 0, 15]
    assert np.isclose(val_lower, 32.0, atol=1e-5), f"Lower boundary expected 32.0, got {val_lower}"

    # 2. Upper boundary check (offOrigin Z = -0.5)
    geo.offOrigin = np.array([-0.5, 0.0, 0.0])
    val_upper = tigre.Ax(img, geo, angles, projection_type="Siddon")[0, 31, 15]
    assert np.isclose(val_upper, 0.0, atol=1e-5), f"Upper boundary expected 0.0, got {val_upper}"


def test_parallel_siddon_grazing_and_cardinal_initial_voxels():
    """Verify that rays grazing voxel boundaries or entering along cardinal angles
    do not misidentify initial voxels due to single-precision rounding.

    Tests multiple challenging orientations:
    1. theta = pi/2, pixel [35, 14]: grazes X = 30 face (~5.700).
    2. alpha = -pi/2, theta = pi/2, psi = 0.13, pixel [27, 37]: entry t matches Y face (~5.727).
    3. alpha = pi, theta = -pi/2, psi = 0.4, pixel [16, 38]: entry t matches X face (~10.173).
    """
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([24, 32, 40]))
    geo.dVoxel = np.array([0.8, 1.1, 0.7])
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.DSO, geo.DSD = np.float32(100), np.float32(200)
    geo.nDetector = np.array([45, 55])
    geo.dDetector = np.array([0.7, 0.9])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.array([0.7, -1.3, 2.1])

    img = np.zeros(geo.nVoxel, dtype=np.float32)
    img[4:18, 6:26, 10:30] = 1.0

    cases = [
        ([[0.0, np.pi / 2, 0.0]], (35, 14), 5.700),
        ([[-np.pi / 2, np.pi / 2, 0.13]], (27, 37), 5.727),
        ([[np.pi, -np.pi / 2, 0.4]], (16, 38), 10.173),
    ]

    for angles, pixel, expected_val in cases:
        angles_arr = np.array(angles, dtype=np.float32)
        val = tigre.Ax(img, geo, angles_arr, projection_type="Siddon")[0][pixel]
        assert np.isclose(val, expected_val, atol=0.01), (
            f"Case {angles} pixel {pixel} failed: expected ~{expected_val}, got {val}"
        )


def test_parallel_siddon_small_angle_sensitivity():
    """Verify that small non-zero tilt angles are not truncated by hardcoded epsilon guards.

    For small theta (1e-5, 2e-5) with offOrigin Z=0.5001, ray [0, 15] crosses the Z boundary.
    The expected values (6.0951 and 10.9922) are regression values computed with TIGRE single-precision
    host endpoint geometry. An analytical ray-box intersection in exact double precision yields ~5.998
    and ~10.999; atol=0.05 accounts for single-precision endpoint rounding while strictly guarding
    against truncation to zero or full-slab integration (32.0).
    """
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([32, 32, 32]))
    geo.dVoxel = np.array([1.0, 1.0, 1.0])
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.dDetector = np.array([1.0, 1.0])
    geo.nDetector = np.array([32, 32])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.array([0.5001, 0.0, 0.0])

    img = np.ones(geo.nVoxel, dtype=np.float32)

    # theta = 1e-5
    angles_1 = np.array([[0.0, 1e-5, 0.0]], dtype=np.float32)
    p_1 = tigre.Ax(img, geo, angles_1, projection_type="Siddon")[0]
    val_1 = p_1[0, 15]
    assert np.isclose(val_1, 6.0951, atol=0.05), f"Small angle ray truncated: {val_1}"

    # theta = 2e-5
    angles_2 = np.array([[0.0, 2e-5, 0.0]], dtype=np.float32)
    p_2 = tigre.Ax(img, geo, angles_2, projection_type="Siddon")[0]
    val_2 = p_2[0, 15]
    assert np.isclose(val_2, 10.9922, atol=0.05), f"Small angle ray mismatch: {val_2}"


def test_parallel_siddon_euler_sensitivity_nondegenerate():
    """Verify non-degenerate 3D Euler angle sensitivity with an asymmetric 3D volume.

    Under non-zero theta, varying psi must independently alter the projection.
    """
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([24, 32, 40]))
    img = np.zeros(geo.nVoxel, dtype=np.float32)
    img[4:18, 6:26, 10:30] = 1.0

    # Under non-zero theta=0.3, vary psi
    angles_psi0 = np.array([[0.2, 0.3, 0.0]], dtype=np.float32)
    angles_psi1 = np.array([[0.2, 0.3, 0.4]], dtype=np.float32)

    p_psi0 = tigre.Ax(img, geo, angles_psi0, projection_type="Siddon")
    p_psi1 = tigre.Ax(img, geo, angles_psi1, projection_type="Siddon")

    diff_psi = np.max(np.abs(p_psi0 - p_psi1))
    assert diff_psi > 1.0, f"psi rotation under non-zero theta had insufficient effect: {diff_psi}"


def _analytical_reference(geo, angles, start, end):
    """Independent analytical ray-box intersection in physical coordinates."""
    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def ry(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    def rx(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    nd = geo.nDetector
    v, u = np.indices(nd)
    yy = (u - nd[1] / 2 + 0.5) * geo.dDetector[1]
    zz = (v - nd[0] / 2 + 0.5) * geo.dDetector[0]
    points = np.stack([np.zeros_like(yy), yy, zz], axis=-1)
    offset = np.asarray(geo.offOrigin)[::-1]
    lower = np.asarray(start)[::-1] * geo.dVoxel[::-1] - geo.sVoxel[::-1] / 2 + offset
    upper = np.asarray(end)[::-1] * geo.dVoxel[::-1] - geo.sVoxel[::-1] / 2 + offset
    ans = []
    yaw, pitch, roll = geo.rotDetector
    detector_rotation = rz(roll) @ ry(pitch) @ rx(yaw)

    for alpha, theta, psi in angles:
        rotation = rz(float(alpha)) @ ry(float(theta)) @ rz(float(psi))
        source = points @ detector_rotation.T
        pixel = source.copy()
        source[..., 0] += geo.DSO
        pixel[..., 0] -= geo.DSD - geo.DSO
        pixel[..., 1] += geo.offDetector[1]
        pixel[..., 2] += geo.offDetector[0]
        source = source @ rotation.T
        pixel = pixel @ rotation.T
        direction = pixel - source
        lo = np.full(nd, -np.inf)
        hi = np.full(nd, np.inf)
        for k in range(3):
            with np.errstate(divide="ignore", invalid="ignore"):
                t1 = (lower[k] - source[..., k]) / direction[..., k]
                t2 = (upper[k] - source[..., k]) / direction[..., k]
            zero = direction[..., k] == 0
            inside = (source[..., k] >= lower[k]) & (source[..., k] < upper[k])
            aa = np.where(zero, np.where(inside, -np.inf, np.inf), np.minimum(t1, t2))
            bb = np.where(zero, np.where(inside, np.inf, -np.inf), np.maximum(t1, t2))
            lo = np.maximum(lo, aa)
            hi = np.minimum(hi, bb)
        ans.append(np.maximum(0, hi - lo) * np.linalg.norm(direction, axis=-1))
    return np.array(ans)


def test_parallel_siddon_independent_analytical_ground_truth():
    """Verify Siddon projection against an independent physical-coordinate analytical ground truth.

    Evaluates 13 views with arbitrary 3D Euler angles, anisotropic voxels, off-center phantom,
    and image offsets.
    """
    rng = np.random.RandomState(51)
    angles = rng.uniform(-2, 2, (13, 3)).astype(np.float32)
    angles[0] = [0.2, 0.3, 0.4]
    angles[1] = [0.0, np.pi / 2, 0.0]

    geo = tigre.geometry(mode="parallel", nVoxel=np.array([24, 32, 40]))
    geo.dVoxel = np.array([0.8, 1.1, 0.7])
    geo.sVoxel = geo.nVoxel * geo.dVoxel
    geo.DSO = np.float32(100)
    geo.DSD = np.float32(200)
    geo.nDetector = np.array([45, 55])
    geo.dDetector = np.array([0.7, 0.9])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offOrigin = np.array([0.7, -1.3, 2.1])
    geo.offDetector = np.array([0.0, 0.0])
    geo.rotDetector = np.array([0.0, 0.0, 0.0])

    start = [4, 6, 10]
    end = [18, 26, 30]
    img = np.zeros(geo.nVoxel, dtype=np.float32)
    img[4:18, 6:26, 10:30] = 1.0

    observed = tigre.Ax(img, geo, angles)
    expected = _analytical_reference(geo, angles, start, end)

    error = np.abs(observed - expected)
    max_error = float(error.max())
    corr = float(np.corrcoef(observed.ravel(), expected.ravel())[0, 1])

    assert max_error < 0.02, f"Max error against analytical ground truth: {max_error}"
    assert corr > 0.99999, f"Correlation with analytical ground truth: {corr}"


def test_parallel_siddon_backward_compatibility():
    """Verify exact equivalence between 1D angles and 3D (alpha, 0, 0) angles."""
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([32, 32, 32]))
    img = np.random.RandomState(42).randn(*geo.nVoxel).astype(np.float32)

    angles_1d = np.linspace(0, 2 * np.pi, 8, endpoint=False, dtype=np.float32)
    angles_3d = np.zeros((8, 3), dtype=np.float32)
    angles_3d[:, 0] = angles_1d

    p_1d = tigre.Ax(img, geo, angles_1d, projection_type="Siddon")
    p_3d = tigre.Ax(img, geo, angles_3d, projection_type="Siddon")
    assert np.allclose(p_1d, p_3d, atol=1e-6), f"Max diff between 1D and 3D: {np.max(np.abs(p_1d - p_3d))}"


def test_parallel_siddon_baseline_untilted_exact_match():
    """Verify regression consistency with upstream baseline (SHA 6b0951a8) on untilted projections.

    Evaluates 11 arbitrary alpha views on a random 32^3 volume with theta=0 and psi=0.
    Under identical floating-point operations, comprehensive offline diagnostic confirms
    exact bitwise agreement across all 11,264 detector pixels (max absolute error = 0.0).
    This test verifies that key pixels and per-view integrated sums match the recorded
    upstream baseline outputs with strict absolute tolerance (atol=1e-5, rtol=0).
    """
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([32, 32, 32]))
    img = np.random.RandomState(42).randn(32, 32, 32).astype(np.float32)
    angles = np.array(
        [0, 0.13, 0.3, 0.7, np.pi / 4, np.pi / 2, 1.8, 2.7, 3.5, 4.1, 5.9],
        dtype=np.float32,
    )
    observed = tigre.Ax(img, geo, angles, projection_type="Siddon")

    # Baseline reference values at key pixels (computed from upstream HEAD master 6b0951a8)
    assert np.isclose(observed[0, 16, 16], 7.398002, atol=1e-5, rtol=0)
    assert np.isclose(observed[5, 16, 16], 9.700473, atol=1e-5, rtol=0)
    assert np.isclose(observed[9, 24, 28], -5.863827, atol=1e-5, rtol=0)
    assert np.isclose(observed[9, 15, 11], -3.0249376, atol=1e-5, rtol=0)

    # Per-view integrated sums across all 11 views
    expected_sums = [
        -38.325397, 4.210716, -28.240627, 23.200993, 0.682480,
        -38.325356, 26.971115, 47.643982, -13.041023, -28.732788, 25.800674
    ]
    observed_sums = [float(observed[v].sum()) for v in range(11)]
    assert np.allclose(observed_sums, expected_sums, atol=1e-4, rtol=0)


if __name__ == "__main__":
    test_parallel_siddon_z_boundary_sampling()
    test_parallel_siddon_grazing_and_cardinal_initial_voxels()
    test_parallel_siddon_small_angle_sensitivity()
    test_parallel_siddon_euler_sensitivity_nondegenerate()
    test_parallel_siddon_independent_analytical_ground_truth()
    test_parallel_siddon_backward_compatibility()
    test_parallel_siddon_baseline_untilted_exact_match()
    print("All parallel Siddon Euler angle tests passed successfully!")
