"""Ax and Atb must not change the GpuIds they are given (#855).

When there are fewer angles (Ax) or image slices (Atb) than GPUs, they use
fewer GPUs. They used to do that by trimming the caller's own GpuIds, so an
algorithm that projects one angle at a time (SART, ASD_POCS, ...) ran on one
GPU after its first call, for the projections and for the TV step.

The compiled operators are replaced by fakes that record the GPUs they get,
so these tests need neither a GPU nor more than one.
"""
import importlib

import numpy as np

import tigre
from tigre.utilities.gpu import GpuIds

ax_mod = importlib.import_module("tigre.utilities.Ax")
atb_mod = importlib.import_module("tigre.utilities.Atb")


def _two_gpus():
    return GpuIds._from_devices([0, 1])


def test_ax_trims_a_copy(monkeypatch):
    seen = []

    def fake_ax(img, geo, angles, projection_type, mode, gpuids=None):
        seen.append(list(gpuids.devices))
        return np.zeros((len(angles), *geo.nDetector), dtype=np.float32)

    monkeypatch.setattr(ax_mod, "_Ax_ext", fake_ax)
    geo = tigre.geometry(mode="cone", nVoxel=np.array([8, 8, 8]), default=True)
    img = np.zeros((8, 8, 8), dtype=np.float32)
    ids = _two_gpus()

    ax_mod.Ax(img, geo, np.array([0.0]), gpuids=ids)  # one angle: one GPU
    ax_mod.Ax(img, geo, np.linspace(0, np.pi, 4), gpuids=ids)  # four angles: both

    assert seen == [[0], [0, 1]]
    assert ids.devices == [0, 1]


def test_atb_trims_a_copy(monkeypatch):
    seen = []

    def fake_atb(proj, geo, angles, backprojection_type, mode, gpuids=None):
        seen.append(list(gpuids.devices))
        return np.zeros(geo.nVoxel, dtype=np.float32)

    monkeypatch.setattr(atb_mod, "_Atb_ext", fake_atb)
    geo = tigre.geometry(mode="cone", nVoxel=np.array([1, 8, 8]), default=True)
    angles = np.linspace(0, np.pi, 4)
    proj = np.zeros((4, *geo.nDetector), dtype=np.float32)
    ids = _two_gpus()

    atb_mod.Atb(proj, geo, angles, gpuids=ids)  # one slice: one GPU

    assert seen == [[0]]
    assert ids.devices == [0, 1]
