"""CGLS: the float64 norm accumulator, and the restart that never restarted.

Two independent defects, both invisible at demo sizes and both decisive at
production ones:

1. ``np.linalg.norm`` accumulates the sum of squares in the array's OWN dtype.
   On a projection-sized float32 array the running sum saturates against the
   increments still being added to it and the norm comes back LOW - measured
   26 % low at 720 x 3072^2. Every Krylov step size is a ratio of such norms,
   and the loss-of-orthogonality test compares consecutive ones.

2. The port of MATLAB's CGLS.m dropped the OUTER of its two loops, so the
   ``break`` that MATLAB uses to fall back and restart instead left
   ``run_main_iter`` entirely. Together with a sentinel that reads
   ``re_init_at_iteration + 1 == i`` against an initial 0 - true at i == 1 -
   one lost step at the first iteration ended the reconstruction.

Neither test needs a GPU: the numerics are pure numpy, and the control-flow
test drives run_main_iter with scripted residual norms.
"""
import numpy as np
import pytest

from tigre.algorithms import krylov_subspace_algorithms as K


def test_norm_is_accurate_where_float32_accumulation_saturates():
    """Equal-magnitude terms: the accumulator grows while the terms do not.

    The error accelerates with N rather than tracking sqrt(N) - it is
    progressive saturation, not random walk - so this modest size shows the
    effect while the production sizes in the module docstring show its bite.
    """
    n = 50_000_000
    x = np.full(n, 1e-3, dtype=np.float32)
    exact = np.sqrt(n) * 1e-3

    accurate = float(K.l2norm(x))
    naive = float(np.linalg.norm(x))
    print(f"exact {exact:.6f}  helper {accurate:.6f}  np.linalg.norm {naive:.6f}")

    assert accurate == pytest.approx(exact, rel=1e-5)
    # The helper must return float32 so callers' dtypes are unchanged:
    # these scalars multiply float32 volumes and Ax/Atb reject float64.
    assert K.l2norm(x).dtype == np.float32


def test_norm_matches_numpy_when_accumulation_is_safe():
    rng = np.random.default_rng(0)
    x = rng.standard_normal(10_000, dtype=np.float32)
    assert float(K.l2norm(x)) == pytest.approx(float(np.linalg.norm(x)), rel=1e-6)
    # Non-float32 input is passed straight through.
    y = rng.standard_normal(1000)
    assert float(K.l2norm(y)) == pytest.approx(float(np.linalg.norm(y)), rel=1e-12)


def _stub_cgls(monkeypatch, residuals, niter, restart=True):
    """A CGLS instance whose residual norms follow `residuals`, no GPU.

    Per inner iteration the projection-shaped norms are taken in the order
    (q_norm, residual); volume-shaped norms are p_norm / s_norm. Returning 1.0
    for everything except the residual keeps alpha and beta finite and makes
    the residual sequence the only thing driving the control flow.
    """
    proj_shape, vol_shape = (4, 5, 5), (3, 3, 3)
    alg = object.__new__(K.CGLS)
    alg.niter = niter
    alg.verbose = False
    alg.Quameasopts = None
    alg.restart = restart
    alg.re_init_at_iteration = -1
    alg.geo = None
    alg.angles = np.zeros(proj_shape[0], dtype=np.float32)
    alg.gpuids = None
    alg.proj = np.zeros(proj_shape, dtype=np.float32)
    alg.res = np.zeros(vol_shape, dtype=np.float32)

    seen = {"proj_norms": 0, "restarts": 0, "residuals": []}

    def fake_norm(x, ord=2):
        x = np.asarray(x)
        if x.size == int(np.prod(proj_shape)):
            j = seen["proj_norms"]
            seen["proj_norms"] += 1
            if j % 2 == 0:                      # q_norm
                return np.float32(1.0)
            k = j // 2                          # residual for iteration k
            v = residuals[min(k, len(residuals) - 1)]
            seen["residuals"].append(v)
            return np.float32(v)
        return np.float32(1.0)

    def fake_Ax(x, geo, angles, ktype, gpuids=None):
        return np.zeros(proj_shape, dtype=np.float32)

    def fake_Atb(x, geo, angles, backprojection_type=None, gpuids=None):
        return np.zeros(vol_shape, dtype=np.float32)

    monkeypatch.setattr(K, "l2norm", fake_norm)
    monkeypatch.setattr(K, "Ax", fake_Ax)
    monkeypatch.setattr(K, "Atb", fake_Atb)
    monkeypatch.setattr(K.tigre, "Ax", fake_Ax)
    monkeypatch.setattr(K.tigre, "Atb", fake_Atb)

    real_init = alg.initialize_algo
    def counting_init():
        seen["restarts"] += 1
        return real_init()
    alg.initialize_algo = counting_init
    return alg, seen


def test_a_lost_step_restarts_instead_of_ending_the_reconstruction(monkeypatch):
    """Residual rises once at iteration 1, then falls: must run all niter."""
    residuals = [10.0, 11.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0]
    alg, seen = _stub_cgls(monkeypatch, residuals, niter=6)
    alg.run_main_iter()

    # One rebuild at entry, one more for the restart - and crucially the run
    # did not stop at i == 1: every iteration slot was filled.
    assert seen["restarts"] == 2
    assert np.count_nonzero(alg.l2l[0]) == 6


def test_repeated_failure_at_the_same_iteration_still_gives_up(monkeypatch):
    """The restart is allowed one chance, exactly as in the MATLAB reference."""
    alg, seen = _stub_cgls(monkeypatch, [10.0, 11.0, 11.0, 11.0], niter=6)
    alg.run_main_iter()
    assert seen["restarts"] == 2            # entry + one restart, then exit
    assert np.count_nonzero(alg.l2l[0]) < 6


def test_restart_false_preserves_the_early_return(monkeypatch):
    alg, seen = _stub_cgls(monkeypatch, [10.0, 11.0, 8.0], niter=6, restart=False)
    alg.run_main_iter()
    assert seen["restarts"] == 1
    assert np.count_nonzero(alg.l2l[0]) == 2


def test_im3dnorm_l2_uses_the_same_accumulator():
    """ASD-POCS reaches the accumulator through im3DNORM, not via the Krylov module.

    Its data-fit test and TV step control are ratios of these norms
    (`dd` is projection-sized, `dp`/`dg` volume-sized), so the two entry points
    must not drift apart.
    """
    from tigre.utilities.im3Dnorm import im3DNORM, l2norm
    n = 50_000_000
    x = np.full(n, 1e-3, dtype=np.float32)
    exact = np.sqrt(n) * 1e-3
    assert float(im3DNORM(x, 2)) == pytest.approx(exact, rel=1e-5)
    assert K.l2norm is l2norm
    # Other norms are untouched.
    y = np.arange(10, dtype=np.float32)
    assert float(im3DNORM(y, 1)) == pytest.approx(float(np.linalg.norm(y, 1)))


def test_inner_is_accurate_and_keeps_float32():
    """Gram-Schmidt coefficients go through inner(); their error compounds
    across a Krylov basis, so they get the float64 accumulator too."""
    from tigre.utilities.im3Dnorm import inner
    n = 50_000_000
    a = np.full(n, 1e-3, dtype=np.float32)
    b = np.full(n, 2e-3, dtype=np.float32)
    exact = n * 2e-6
    assert float(inner(a, b)) == pytest.approx(exact, rel=1e-5)
    assert inner(a, b).dtype == np.float32
    # Shapes are ravelled, so a volume and its flat view agree.
    v = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    assert float(inner(v, v)) == pytest.approx(float(np.dot(v.ravel(), v.ravel())), rel=1e-6)


def test_lsqr_and_lsmr_share_the_restart_fix():
    """They had the identical single-loop structure and the identical
    sentinel; a fix that reached only CGLS would leave them stopping at
    iteration 1."""
    import inspect
    for cls_name in ("LSQR", "LSMR", "CGLS"):
        src = inspect.getsource(getattr(K, cls_name).run_main_iter)
        assert src.count("while i < self.niter:") == 2, f"{cls_name}: outer loop missing"
        # The code form, not the docstrings - CGLS's explains the old bug.
        assert "if self.re_init_at_iteration + 1 == i" not in src, f"{cls_name}: old sentinel"
        assert "if self.re_init_at_iteration == i" in src, f"{cls_name}: give-up test missing"


def test_irn_tv_difference_operator_is_a_true_adjoint_pair():
    """IRN_TV_CGLS runs CGLS on the stacked operator [A; sqrt(lmbda) W D],
    which is only meaningful if the transpose it is handed IS the adjoint.

    The original pair was wrong twice over: `Dxx = np.copy(img)` followed by
    writing only `[0:-2]` left the last two slices holding raw image VALUES
    instead of differences - so D did not annihilate a constant - and D^T was
    off by one at index n-2. Measured asymmetry 6.4e-2, and the inner CGLS
    diverged to 7e5 on a 64^3 phantom.
    """
    from tigre.algorithms.krylov_subspace_algorithms import IRN_TV_CGLS

    alg = object.__new__(IRN_TV_CGLS)
    rng = np.random.default_rng(0)
    n = 12
    W = (rng.random((n, n, n)).astype(np.float32) + 0.5)
    u = rng.standard_normal((n, n, n), dtype=np.float32)
    v = rng.standard_normal((3, n, n, n), dtype=np.float32)

    lhs = float(np.sum(IRN_TV_CGLS.Lx(alg, W, u) * v))
    rhs = float(np.sum(u * IRN_TV_CGLS.Ltx(alg, W, v)))
    assert lhs == pytest.approx(rhs, rel=1e-5), f"<Lu,v>={lhs} vs <u,Ltv>={rhs}"

    # A gradient operator annihilates a constant image, at the boundary too.
    ones = np.ones((n, n, n), dtype=np.float32)
    assert np.abs(IRN_TV_CGLS.Lx(alg, ones, ones)).max() == 0.0
