"""End-to-end GPU tests: every reconstruction algorithm on the head phantom.

Each algorithm reconstructs the 64^3 head phantom from noise-free projections and
must return a finite volume of the right shape whose relative L2 error against the
phantom is below a bound. Idea for one test per algorithm from PR #768 (Vilax).

The bounds in MAX_REL_ERROR are deliberately loose to start with; they catch an
algorithm that diverges, returns garbage or stops working, not small regressions.
Tighten them once you have numbers from your own GPU: run with `-s` and each test
prints the error it measured.

Run them locally on a machine with a GPU:  pytest -m gpu Python/tests
"""
import numpy as np
import pytest

import tigre
import tigre.algorithms as algs
from tigre.utilities import sample_loader

pytestmark = pytest.mark.gpu

N_ANGLES = 60

# name, iterations, extra keyword arguments
ITERATIVE = [
    ("sart", 3, {}),
    ("sirt", 30, {}),
    ("ossart", 10, {}),
    ("sart_tv", 3, {}),
    ("ossart_tv", 10, {}),
    ("cgls", 15, {}),
    ("lsqr", 15, {}),
    ("hybrid_lsqr", 15, {}),
    ("lsmr", 15, {}),
    ("irn_tv_cgls", 10, {"lmbda": 5, "niter_outer": 2}),
    ("hybrid_flsqr_tv", 15, {}),
    ("ab_gmres", 15, {}),
    ("ba_gmres", 15, {}),
    ("asd_pocs", 10, {}),
    ("os_asd_pocs", 10, {}),
    ("awasd_pocs", 10, {}),
    ("os_awasd_pocs", 10, {}),
    ("pcsd", 10, {}),
    ("aw_pcsd", 10, {}),
    ("os_pcsd", 10, {}),
    ("os_aw_pcsd", 10, {}),
    ("mlem", 30, {}),
    ("fista", 30, {}),
    ("ista", 30, {}),
]

PRIOR_BASED = [("piccs", 10), ("os_piccs", 10)]

# Relative L2 error bound per algorithm (default for anything not listed).
DEFAULT_MAX_REL_ERROR = 0.6
MAX_REL_ERROR = {}


def rel_error(x, truth):
    return float(np.linalg.norm(x - truth) / np.linalg.norm(truth))


def result_volume(out):
    # Algorithms return (image, errors) when computel2 or Quameasopts is set.
    return out[0] if isinstance(out, tuple) else out


@pytest.fixture(scope="module")
def cone_problem():
    geo = tigre.geometry(mode="cone", default=True, high_resolution=False)
    angles = np.linspace(0, 2 * np.pi, N_ANGLES, endpoint=False, dtype=np.float32)
    truth = sample_loader.load_head_phantom(geo.nVoxel).astype(np.float32)
    proj = tigre.Ax(truth, geo, angles)
    return geo, angles, truth, proj


@pytest.fixture(scope="module")
def parallel_problem():
    geo = tigre.geometry(mode="parallel", nVoxel=np.array([64, 64, 64]))
    angles = np.linspace(0, np.pi, N_ANGLES, endpoint=False, dtype=np.float32)
    truth = sample_loader.load_head_phantom(geo.nVoxel).astype(np.float32)
    proj = tigre.Ax(truth, geo, angles)
    return geo, angles, truth, proj


def check(name, img, geo, truth):
    assert img.shape == tuple(geo.nVoxel), name
    assert np.isfinite(img).all(), name
    err = rel_error(img, truth)
    print(f"\n{name}: relative L2 error {err:.4f}")
    assert err < MAX_REL_ERROR.get(name, DEFAULT_MAX_REL_ERROR), f"{name}: error {err:.4f}"


def test_fdk(cone_problem):
    geo, angles, truth, proj = cone_problem
    check("fdk", algs.fdk(proj, geo, angles), geo, truth)


@pytest.mark.parametrize("filt", ["ram_lak", "shepp_logan", "cosine", "hamming", "hann"])
def test_fdk_filters(cone_problem, filt):
    geo, angles, truth, proj = cone_problem
    check(f"fdk[{filt}]", algs.fdk(proj, geo, angles, filter=filt), geo, truth)


def test_fbp_parallel(parallel_problem):
    geo, angles, truth, proj = parallel_problem
    check("fbp", algs.fbp(proj, geo, angles), geo, truth)


@pytest.mark.parametrize("name,niter,kwargs", ITERATIVE, ids=[a[0] for a in ITERATIVE])
def test_iterative_cone(cone_problem, name, niter, kwargs):
    geo, angles, truth, proj = cone_problem
    out = getattr(algs, name)(proj, geo, angles, niter, verbose=False, **kwargs)
    check(name, result_volume(out), geo, truth)


@pytest.mark.parametrize("name,niter", PRIOR_BASED, ids=[a[0] for a in PRIOR_BASED])
def test_prior_based_cone(cone_problem, name, niter):
    geo, angles, truth, proj = cone_problem
    prior = algs.fdk(proj, geo, angles)
    out = getattr(algs, name)(proj, geo, angles, niter, prior=prior, verbose=False)
    check(name, result_volume(out), geo, truth)


@pytest.mark.parametrize("name,niter", [("sirt", 30), ("ossart", 10), ("cgls", 15)])
def test_iterative_parallel(parallel_problem, name, niter):
    geo, angles, truth, proj = parallel_problem
    out = getattr(algs, name)(proj, geo, angles, niter, verbose=False)
    check(f"{name}[parallel]", result_volume(out), geo, truth)


@pytest.mark.parametrize("name", ["cgls", "sirt"])
def test_more_iterations_fit_the_data_better(cone_problem, name):
    """CGLS and SIRT reduce the data residual ||Ax - b|| as iterations go on."""
    geo, angles, truth, proj = cone_problem

    def residual(niter):
        img = result_volume(getattr(algs, name)(proj, geo, angles, niter, verbose=False))
        return np.linalg.norm(tigre.Ax(img, geo, angles) - proj)

    assert residual(20) < residual(5)
