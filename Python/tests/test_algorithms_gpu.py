"""End-to-end GPU tests: every reconstruction algorithm on the head phantom.

Each algorithm reconstructs the 64^3 head phantom from noise-free projections and
must return a finite volume of the right shape whose relative L2 error against the
phantom is below a bound. Idea for one test per algorithm from PR #768 (Vilax).

The bounds in MAX_REL_ERROR are the errors measured on a first run plus some
headroom, so they catch an algorithm that gets noticeably worse. Run with `-s` to
see the error each test measures.

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

# Relative L2 error bound per algorithm: the error measured on 2026-11-11
# (Windows, CUDA 12.3, 1 GPU) plus 25% and 0.01 of headroom. A failure here means
# an algorithm got noticeably worse, not just different rounding. Anything not
# listed falls back to DEFAULT_MAX_REL_ERROR.
DEFAULT_MAX_REL_ERROR = 0.6
MAX_REL_ERROR = {
    "fdk": 0.21,  # measured 0.1521
    "fdk[ram_lak]": 0.21,  # measured 0.1521
    "fdk[shepp_logan]": 0.19,  # measured 0.1387
    "fdk[cosine]": 0.18,  # measured 0.1303
    "fdk[hamming]": 0.18,  # measured 0.1335
    "fdk[hann]": 0.19,  # measured 0.1362
    "fbp": 0.17,  # measured 0.1228
    "sart": 0.09,  # measured 0.0640
    "sirt": 0.28,  # measured 0.2131
    "ossart": 0.28,  # measured 0.2126
    "sart_tv": 0.13,  # measured 0.0894
    "ossart_tv": 0.33,  # measured 0.2549
    "cgls": 0.11,  # measured 0.0775
    "lsqr": 0.12,  # measured 0.0805
    "hybrid_lsqr": 0.12,  # measured 0.0805
    "lsmr": 0.13,  # measured 0.0886
    "irn_tv_cgls": 0.08,  # measured 0.0556
    # GCV (as in MATLAB) over-regularises this noise-free problem: lambda ~120
    # in Python, ~250 in MATLAB's own hybrid_fLSQR_TV (error 0.37 there).
    # With a fixed lmbda=1 the error is 0.177.
    "hybrid_flsqr_tv": 0.45,  # measured 0.3378
    "ab_gmres": 0.12,  # measured 0.0802
    "ba_gmres": 0.12,  # measured 0.0880
    "asd_pocs": 0.06,  # measured 0.0364
    "os_asd_pocs": 0.31,  # measured 0.2356
    "awasd_pocs": 0.04,  # measured 0.0164
    "os_awasd_pocs": 0.29,  # measured 0.2196
    "pcsd": 0.04,  # measured 0.0227
    "aw_pcsd": 0.03,  # measured 0.0136
    "os_pcsd": 0.31,  # measured 0.2383
    "os_aw_pcsd": 0.28,  # measured 0.2134
    "mlem": 0.13,  # measured 0.0925
    "fista": 0.11,  # measured 0.0727
    "ista": 0.20,  # measured 0.1443
    "piccs": 0.05,  # measured 0.0262
    "os_piccs": 0.29,  # measured 0.2200
    "sirt[parallel]": 0.23,  # measured 0.1700
    "ossart[parallel]": 0.22,  # measured 0.1677
    "cgls[parallel]": 0.13,  # measured 0.0901
}


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
