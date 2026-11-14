"""Algorithm logic on the CPU, with Ax/Atb replaced by a small dense matrix.

TIGRE's projectors need a GPU, but the algorithms around them are plain NumPy.
Replacing the projectors by an explicit matrix A (and its exact transpose) lets
these tests compare the algorithms with direct linear-algebra solutions.
"""
import numpy as np
import pytest

import tigre
import tigre.algorithms as algs
from tigre.algorithms import iterative_recon_alg as base_mod
from tigre.algorithms import krylov_subspace_algorithms as krylov_mod

N_VOXEL = np.array([2, 3, 3])
N_DET = np.array([2, 3])
N_ANGLES = 5


@pytest.fixture
def problem(monkeypatch):
    rng = np.random.default_rng(1)
    n = int(np.prod(N_VOXEL))
    m = N_ANGLES * int(np.prod(N_DET))
    A = rng.random((m, n)) + np.eye(m, n)  # well conditioned, m > n

    def Ax(img, geo, angles, projection_type="Siddon", **kwargs):
        out = A @ np.asarray(img, dtype=np.float64).ravel()
        return out.reshape(len(angles), *N_DET).astype(np.float32)

    def Atb(proj, geo, angles, backprojection_type="FDK", **kwargs):
        out = A.T @ np.asarray(proj, dtype=np.float64).ravel()
        return out.reshape(N_VOXEL).astype(np.float32)

    for mod in (base_mod, krylov_mod):
        monkeypatch.setattr(mod, "Ax", Ax)
        monkeypatch.setattr(mod, "Atb", Atb)
    monkeypatch.setattr(tigre, "Ax", Ax)
    monkeypatch.setattr(tigre, "Atb", Atb)

    geo = tigre.geometry(mode="cone", nVoxel=N_VOXEL.copy(), default=True)
    geo.nDetector = N_DET.copy()
    geo.dDetector = np.array([1.0, 1.0])
    geo.sDetector = geo.nDetector * geo.dDetector
    angles = np.linspace(0, 2 * np.pi, N_ANGLES, endpoint=False).astype(np.float32)
    x_true = rng.random(n)
    b = A @ x_true
    proj = b.reshape(N_ANGLES, *N_DET).astype(np.float32)
    return A, geo, angles, proj, x_true


def tikhonov(A, b, lmbda):
    n = A.shape[1]
    return np.linalg.solve(A.T @ A + lmbda**2 * np.eye(n), A.T @ b)


def test_lsmr_default_is_undamped_least_squares(problem):
    A, geo, angles, proj, x_true = problem
    x = algs.lsmr(proj, geo, angles, A.shape[1], verbose=False)
    np.testing.assert_allclose(x.ravel(), x_true, rtol=1e-3, atol=1e-3)


def test_lsmr_lambda_sets_damping(problem):
    A, geo, angles, proj, _ = problem
    alg = krylov_mod.LSMR(proj, geo, angles, 2, lmbda=3.0, verbose=False)
    assert alg.lmbda == 3.0
    assert krylov_mod.LSMR(proj, geo, angles, 2, verbose=False).lmbda == 0


def test_hybrid_lsqr_given_lambda_is_tikhonov(problem):
    A, geo, angles, proj, _ = problem
    x = algs.hybrid_lsqr(proj, geo, angles, A.shape[1], lmbda=2.0, verbose=False)
    np.testing.assert_allclose(x.ravel(), tikhonov(A, proj.ravel().astype(float), 2.0), rtol=1e-3, atol=1e-3)


def test_gcv_choice_matches_brute_force():
    rng = np.random.default_rng(2)
    B = rng.random((9, 8)) + np.eye(9, 8)
    b = rng.random(9)
    U, s, _ = np.linalg.svd(B)
    bhat = U.T @ b
    lam = krylov_mod._choose_lambda("gcv", None, None, B, b, bhat, s, np.linalg.norm(b))
    grid = np.linspace(1e-6, s[0], 200001)
    best = grid[np.argmin([krylov_mod._tikhonov_gcv(g, bhat, s) for g in grid])]
    assert lam == pytest.approx(best, rel=1e-2, abs=1e-4)


def test_hybrid_lsqr_gcv_runs(problem):
    A, geo, angles, proj, _ = problem
    rng = np.random.default_rng(2)
    noisy = (proj.ravel() + 0.05 * rng.standard_normal(proj.size)).reshape(proj.shape).astype(np.float32)
    alg = krylov_mod.hybrid_LSQR(noisy, geo, angles, A.shape[1], verbose=False)
    alg.run_main_iter()
    done = alg.lambda_vec[alg.l2l[0] > 0]
    assert len(done) >= 3 and np.all(done > 0)
    rel_res = np.linalg.norm(A @ alg.getres().ravel() - noisy.ravel()) / np.linalg.norm(noisy)
    assert rel_res < 0.1


def test_hybrid_lsqr_noise_level_meets_discrepancy(problem):
    A, geo, angles, proj, _ = problem
    rng = np.random.default_rng(3)
    b = proj.ravel().astype(float)
    noise = rng.standard_normal(b.size)
    noise *= 0.05 * np.linalg.norm(b) / np.linalg.norm(noise)
    noisy = (b + noise).reshape(proj.shape).astype(np.float32)
    x = algs.hybrid_lsqr(noisy, geo, angles, A.shape[1], noiselevel=0.05, verbose=False)
    rel_res = np.linalg.norm(A @ x.ravel() - noisy.ravel()) / np.linalg.norm(noisy.ravel())
    assert rel_res == pytest.approx(0.05, rel=0.05)


def test_flsqr_tv_derivative_pair_is_adjoint():
    rng = np.random.default_rng(4)
    u = rng.random((3, 4, 5))
    y = rng.random((3, 4, 5))
    for ax in range(3):
        lhs = np.sum(krylov_mod._forward_difference(u, ax) * y)
        rhs = np.sum(u * krylov_mod._forward_difference_adjoint(y, ax))
        assert lhs == pytest.approx(rhs, rel=1e-12)
        # D annihilates constants
        assert np.all(krylov_mod._forward_difference(np.ones_like(u), ax) == 0)


@pytest.mark.parametrize("kwargs", [{}, {"lmbda": 0.5}, {"noiselevel": 0.05}])
def test_hybrid_flsqr_tv_runs_and_fits_data(problem, kwargs):
    A, geo, angles, proj, x_true = problem
    x = algs.hybrid_flsqr_tv(proj, geo, angles, 6, verbose=False, **kwargs)
    rel_res = np.linalg.norm(A @ x.ravel() - proj.ravel()) / np.linalg.norm(proj.ravel())
    assert np.all(np.isfinite(x))
    assert rel_res < 0.2


def test_irn_tv_cgls_default_lambda_is_ten(problem):
    A, geo, angles, proj, _ = problem
    alg = krylov_mod.IRN_TV_CGLS(proj, geo, angles, 2, verbose=False)
    assert alg.lmbda == 10


def test_sart_family_reduces_lambda(problem, monkeypatch):
    from tigre.algorithms.art_family_algorithms import OS_SART

    # only the relaxation schedule is under test, not the update itself
    monkeypatch.setattr(OS_SART, "art_data_minimizing", lambda self: None)

    A, geo, angles, proj, _ = problem
    # W and V are given, so the setup does not call the projectors with other geometries
    nblocks = int(np.ceil(N_ANGLES / 2))
    alg = OS_SART(proj, geo, angles, 4, blocksize=2, lmbda=1.0, lmbda_red=0.5, verbose=False,
                  W=np.ones(proj.shape, np.float32), V=np.ones((nblocks, *N_VOXEL[1:]), np.float32))
    alg.run_main_iter()
    assert alg.lmbda == pytest.approx(1.0 * 0.5**4)
