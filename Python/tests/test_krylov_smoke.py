"""Every Krylov algorithm reconstructs a small phantom without falling over.

The unit tests beside this file check the arithmetic and the restart control
flow in isolation. This is the coarse net: run each algorithm end to end and
require that it finishes, returns finite values, and correlates with the truth
about as well as FDK does on the same data.

It exists because three defects fixed on 2026-08-24 were all invisible to
unit tests and all produced plausible-looking output:

  * CGLS / LSQR / LSMR exited after a single iteration on the first residual
    rise, because the port dropped MATLAB's outer restart loop and the
    give-up sentinel fired at i == 1;
  * AB/BA-GMRES, hybrid_LSQR and hybrid_fLSQR_TV returned their result by
    `return <expr>` - but decorator() discards that and hands back
    `self.res`, so any early exit silently returned the INITIAL volume
    (zeros, or with the default FDK warm start an unimproved FDK image);
  * hybrid_fLSQR_TV never assigned self.res on ANY path, so its output was
    discarded outright;
  * IRN_TV_CGLS's difference operator and its "transpose" were not an adjoint
    pair (6.4e-2 asymmetry) and its inner CGLS had no divergence guard, so it
    ran to 7e5 on this very phantom.

Needs a GPU. Small on purpose: 64^3 at 50 views, a few seconds per algorithm.
"""
import numpy as np
import pytest

import tigre
from tigre.algorithms import (cgls, lsqr, lsmr, ab_gmres, ba_gmres,
                              hybrid_lsqr, irn_tv_cgls, fdk)
from tigre.algorithms.krylov_subspace_algorithms import CGLS

NITER = 5


@pytest.fixture(scope="module")
def problem():
    geo = tigre.geometry(mode="cone", nVoxel=np.array([64, 64, 64]), default=True)
    angles = np.linspace(0, 2 * np.pi, 50, endpoint=False, dtype=np.float32)

    z, y, x = np.mgrid[:64, :64, :64].astype(np.float32)
    phantom = np.zeros((64, 64, 64), dtype=np.float32)
    phantom[((x - 30) ** 2 + (y - 32) ** 2 + (z - 32) ** 2) < 14 ** 2] = 1.0
    phantom[((x - 40) ** 2 + (y - 38) ** 2 + (z - 30) ** 2) < 6 ** 2] = 2.0

    proj = tigre.Ax(phantom, geo, angles)
    ref = float(np.corrcoef(fdk(proj.copy(), geo, angles).ravel(), phantom.ravel())[0, 1])
    return geo, angles, phantom, proj, ref


@pytest.mark.parametrize("alg", [cgls, lsqr, lsmr, ab_gmres, ba_gmres, hybrid_lsqr,
                                 irn_tv_cgls])
def test_algorithm_runs_and_reconstructs(problem, alg):
    geo, angles, phantom, proj, ref = problem
    res = alg(proj.copy(), geo, angles, niter=NITER, verbose=False)

    assert res.shape == phantom.shape
    assert np.isfinite(res).all(), f"{alg.__name__} produced non-finite voxels"
    assert np.any(res != 0), f"{alg.__name__} returned an all-zero volume"

    r = float(np.corrcoef(res.ravel(), phantom.ravel())[0, 1])
    # Judged against FDK on the same data rather than an absolute number: this
    # geometry/phantom pairing is deliberately crude, and FDK itself only
    # reaches ~0.42 on it. A collapsed reconstruction lands near zero.
    assert r > 0.6 * ref, f"{alg.__name__} correlates {r:.3f} vs FDK's {ref:.3f}"


def test_cgls_does_not_stop_after_one_iteration(problem):
    """The regression itself: CGLS used to give up at i == 1 and report
    divergence, leaving 2 of 20 residual slots filled."""
    geo, angles, phantom, proj, _ = problem
    geo.check_geo(angles)
    alg = CGLS(proj.copy(), geo, angles, niter=20, verbose=False)
    alg.run_main_iter()
    assert np.count_nonzero(alg.l2l) > 2, (
        f"CGLS filled only {np.count_nonzero(alg.l2l)} residual slots of 20")


def test_no_builtin_sum_over_arrays():
    """AB/BA-GMRES orthogonalised with Python's builtin sum() over a CT-sized
    array - element by element, 40x slower than np.dot and less accurate. It
    is the kind of line that reads as harmless, so guard it."""
    import inspect
    from tigre.algorithms import krylov_subspace_algorithms as K
    src = inspect.getsource(K)
    offenders = [ln.strip() for ln in src.splitlines()
                 if "=sum(" in ln.replace(" ", "") and "avgtime" not in ln]
    assert not offenders, offenders
