"""gpuids=None must mean "every visible GPU", not "no GPU".

Every compiled entry point declares `gpuids=None` in its own signature, and
they all convert it through one helper, `convert_to_c_gpuids` in
_gpuUtils.pxd. That helper used to turn None into a VALID pointer holding
m_iCount = 0 and a NULL device list, and hand that to CUDA.

The algorithm classes hide it - they all do
`if self.gpuids is None: self.gpuids = GpuIds()` before calling down - so
tigre.Ax/Atb and every algs.* function looked fine. The kernels reachable
directly did not: `minTV(img, alpha, iters, None)` ran with zero devices and
killed the interpreter outright. No exception, no message, no traceback - the
process simply ended, which is why it reads as "TIGRE fails" rather than as a
bad argument.

NOTE ON FAILURE MODE: a regression here does not fail these tests, it ABORTS
the interpreter running them. pytest reports that as a crashed worker rather
than an assertion - loud, but not a normal failure. That is the nature of the
bug being guarded.
"""
import numpy as np
import pytest

import tigre
from tigre.utilities.Ax import Ax
from tigre.utilities.Atb import Atb


@pytest.fixture(scope="module")
def small():
    geo = tigre.geometry(mode="cone", nVoxel=np.array([32, 32, 32]), default=True)
    angles = np.linspace(0, 2 * np.pi, 30, endpoint=False, dtype=np.float32)
    geo.check_geo(angles)
    img = np.random.default_rng(0).random((32, 32, 32)).astype(np.float32)
    return geo, angles, img, Ax(img, geo, angles, "Siddon")


def test_ax_atb_accept_none(small):
    geo, angles, img, proj = small
    assert np.isfinite(Ax(img, geo, angles, "Siddon", gpuids=None)).all()
    assert np.isfinite(Atb(proj, geo, angles, backprojection_type="matched",
                           gpuids=None)).all()


def test_min_tv_accepts_none(small):
    """The case that killed the interpreter."""
    from _minTV import minTV
    geo, angles, img, proj = small
    out = minTV(img, 1.0, 2, None)
    assert out.shape == img.shape
    assert np.isfinite(out).all()
    assert np.any(out != 0), "denoiser returned an all-zero volume"


def test_awmin_tv_accepts_none(small):
    from _AwminTV import AwminTV
    geo, angles, img, proj = small
    out = AwminTV(img, 1.0, 2, -0.005, None)
    assert out.shape == img.shape
    assert np.isfinite(out).all()
    assert np.any(out != 0)


def test_tv_proximal_and_poisson_accept_none(small):
    """_tvdenoising declares the same default but is not among the built
    extensions; _tv_proximal is the one that ships."""
    from _tv_proximal import tvdenoise
    from _RandomNumberGenerator import add_poisson
    geo, angles, img, proj = small
    out = tvdenoise(img, 2, 1.0, gpuids=None)
    assert np.isfinite(out).all()
    noisy = add_poisson(np.abs(proj) + 1.0, gpuids=None)
    assert np.isfinite(noisy).all()


def test_none_resolves_to_all_visible_gpus():
    """The default is 'all of them', matching what the algorithm classes do
    by hand, not 'the first one' and not 'none'."""
    from tigre.utilities.gpu import GpuIds
    assert len(GpuIds()) >= 1
