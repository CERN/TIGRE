"""Parallel-beam Ax and Atb on several GPUs (#286).

They used to run on the first GPU only, whatever ``gpuids`` said. Now Ax
shares the angles between the GPUs, and Atb shares the image slices, as the
cone-beam operators do. The result must be the one GPU gives: the same numbers
for Ax (each projection is computed exactly as on one GPU), and the same to
rounding for Atb (each GPU's slab has its own origin).

Needs 2 GPUs; skipped otherwise.
"""
import numpy as np
import pytest

pytestmark = pytest.mark.gpu

import tigre  # noqa: E402
from tigre.utilities.gpu import GpuIds  # noqa: E402
from device_kinds import same_kind  # noqa: E402  (a GPU and the CPU device: list the GPU twice)


def _gpus(devices):
    ids = GpuIds()
    ids.devices = list(devices)
    return ids


try:
    _n_gpus = len(GpuIds().devices)
except Exception:  # no GPU or no driver
    _n_gpus = 0
if _n_gpus < 2:
    pytest.skip("needs 2 GPUs", allow_module_level=True)


def _problem(nangles, nvoxel):
    geo = tigre.geometry(mode="parallel", default=True, nVoxel=np.array(nvoxel))
    angles = np.linspace(0, 2 * np.pi, nangles, endpoint=False)
    img = np.random.default_rng(0).random(tuple(nvoxel), dtype=np.float32)
    return geo, angles, img


# 37 angles: not a multiple of the projectors' batch sizes; 2: one angle per GPU.
@pytest.mark.parametrize("nangles", [37, 2])
@pytest.mark.parametrize("kind", ["Siddon", "interpolated"])
def test_ax_two_gpus_match_one(kind, nangles):
    geo, angles, img = _problem(nangles, [65, 64, 72])
    one = tigre.Ax(img, geo, angles, kind, gpuids=_gpus([0]))
    two = tigre.Ax(img, geo, angles, kind, gpuids=_gpus(same_kind([0, 1])))
    np.testing.assert_array_equal(two, one)


# 65 slices: the two GPUs' slabs differ in size; 2: one slice per GPU.
@pytest.mark.parametrize("nvoxel", [[65, 64, 72], [2, 64, 72]])
def test_atb_two_gpus_match_one(nvoxel):
    geo, angles, img = _problem(37, nvoxel)
    proj = tigre.Ax(img, geo, angles, gpuids=_gpus([0]))
    one = tigre.Atb(proj, geo, angles, gpuids=_gpus([0]))
    two = tigre.Atb(proj, geo, angles, gpuids=_gpus(same_kind([0, 1])))
    err = float(np.abs(two - one).max() / np.abs(one).max())
    assert err < 1e-5, err
