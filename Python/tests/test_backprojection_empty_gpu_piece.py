"""Backprojection on several GPUs when the last GPU gets no slices.

The cone backprojectors give each GPU a slab of ceil(nz/G) slices; then the last GPUs can be left with fewer, or
none: 5 slices on 4 GPUs are slabs of 2, 2, 1 and 0. The empty slab's size was computed in unsigned arithmetic and
came out as -1 slices: its launch failed ("invalid argument") and its copy back had a huge size, an illegal memory
access that broke every later CUDA call in the process. Needs 4+ GPUs (with 2 or 3 no slab is ever empty).
"""
import numpy as np
import pytest

import tigre
from tigre.utilities.gpu import GpuIds

pytestmark = pytest.mark.gpu


def _gpus():
    try:
        return len(GpuIds())
    except Exception:
        return 0


def _sizes_with_an_empty_slab(g):
    """Slice counts nz >= g (fewer slices than GPUs use fewer GPUs) where ceil(nz/g) slices per GPU leave the last
    GPU none."""
    return [nz for nz in range(g, 4 * g) if -(-nz // g) * (g - 1) >= nz][:3]


@pytest.mark.skipif(_gpus() < 4, reason="needs 4 or more GPUs")
@pytest.mark.parametrize("btype", ["FDK", "matched"])
def test_all_gpus_match_one_when_a_slab_is_empty(btype):
    g = _gpus()
    sizes = _sizes_with_an_empty_slab(g)
    assert sizes, "no such size for {} GPUs".format(g)
    for nz in sizes:
        geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([nz, 48, 48]))
        angles = np.linspace(0, 2 * np.pi, 40, endpoint=False)
        proj = np.random.default_rng(nz).random((len(angles), *geo.nDetector), dtype=np.float32)
        one = tigre.Atb(proj, geo, angles, btype, gpuids=GpuIds()[0])
        every = tigre.Atb(proj, geo, angles, btype, gpuids=GpuIds())
        np.testing.assert_allclose(every, one, rtol=1e-5, atol=1e-5 * float(np.abs(one).max()),
                                   err_msg="{} slices on {} GPUs".format(nz, g))
