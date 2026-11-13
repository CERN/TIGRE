"""minTV / AwminTV must not depend on what ran on the GPU before them.

TIGRE is compiled with --default-stream=per-thread, so a plain cudaMemset runs
on the per-thread default stream and is not ordered with the streams the CUDA
code creates itself. In GD_TV.cu / GD_AwTV.cu the image buffer was zeroed with
cudaMemset and then filled with cudaMemcpyAsync on another stream; sometimes the
zeroing landed last and the TV step ran on an all-zero image, returning zeros.
Whether it happened depended on timing, i.e. on earlier GPU work.

These tests interleave projections and TV calls many times and require every
result to be non-zero and identical to the first one.
"""
import numpy as np
import pytest

import tigre
from tigre.utilities import sample_loader
from tigre.utilities.gpu import GpuIds

pytestmark = pytest.mark.gpu

REPEATS = 30


@pytest.fixture(scope="module")
def problem():
    geo = tigre.geometry(mode="cone", default=True, high_resolution=False)
    angles = np.linspace(0, 2 * np.pi, 60, endpoint=False, dtype=np.float32)
    head = sample_loader.load_head_phantom(geo.nVoxel).astype(np.float32)
    img = np.random.default_rng(0).random((32, 32, 32)).astype(np.float32)
    return geo, angles, head, img


def _repeat_with_gpu_work_in_between(problem, tv_call):
    geo, angles, head, img = problem
    first = tv_call(img.copy())
    assert np.any(first != 0), "TV returned an all-zero volume on the first call"
    for i in range(REPEATS):
        proj = tigre.Ax(head, geo, angles)
        tigre.Atb(proj, geo, angles)
        out = tv_call(img.copy())
        assert np.any(out != 0), f"TV returned an all-zero volume on repeat {i}"
        np.testing.assert_array_equal(out, first, err_msg=f"TV result changed on repeat {i}")


def test_min_tv_is_independent_of_earlier_gpu_work(problem):
    from _minTV import minTV

    _repeat_with_gpu_work_in_between(problem, lambda x: minTV(x, 1.0, 2, GpuIds()))


def test_awmin_tv_is_independent_of_earlier_gpu_work(problem):
    from _AwminTV import AwminTV

    _repeat_with_gpu_work_in_between(problem, lambda x: AwminTV(x, 1.0, 2, -0.005, GpuIds()))


@pytest.mark.parametrize("btype", ["FDK", "matched"])
def test_backprojection_is_independent_of_earlier_gpu_work(problem, btype):
    geo, angles, head, img = problem
    proj = tigre.Ax(head, geo, angles)
    first = tigre.Atb(proj, geo, angles, btype)
    for i in range(REPEATS):
        tigre.Ax(head, geo, angles, "interpolated")
        np.testing.assert_array_equal(tigre.Atb(proj, geo, angles, btype), first,
                                      err_msg=f"backprojection changed on repeat {i}")
