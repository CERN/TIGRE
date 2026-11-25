"""SYCL build: Ax with little device memory.

TIGRE_SYCL_MEMORY_LIMIT_MB caps the memory the SYCL build plans with, so the
image-splitting path of 1a36447 (Ax) runs on any SYCL device, including the
CPU device in CI. Each case runs in a subprocess, since the limit is read from
the environment.

What a split result should be:
- Ax (interpolated): the unsplit result.
- Ax (Siddon): the unsplit result up to rounding, since the pieces add up in a
  different order (as in test_siddon_split.py for the CUDA build, #859).
"""

import os
import subprocess
import sys

import numpy as np
import pytest

pytestmark = pytest.mark.gpu

_gpuUtils = pytest.importorskip("tigre.utilities.cuda_interface._gpuUtils")
if getattr(_gpuUtils, "getBackend", lambda: "cuda")() != "sycl":
    pytest.skip("SYCL build only (TIGRE_SYCL_MEMORY_LIMIT_MB)", allow_module_level=True)

SCRIPT = r"""
import sys
import numpy as np
import tigre

img = np.random.default_rng(0).random((64, 256, 256), dtype=np.float32)  # 16 MB, 256 kB per slice
geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array(img.shape))
geo.sVoxel = np.array([64.0, 256.0, 256.0])
geo.dVoxel = geo.sVoxel / geo.nVoxel
geo.nDetector = np.array([61, 67])
geo.dDetector = geo.sDetector / geo.nDetector
angles = np.linspace(0, 2 * np.pi, 8, endpoint=False, dtype=np.float32)
calls = {
    "Ax": lambda: tigre.Ax(img, geo, angles, "interpolated"),
    "AxSiddon": lambda: tigre.Ax(img, geo, angles, "Siddon"),
}
out = {}
for name, call in calls.items():
    try:
        out[name] = call()
    except RuntimeError as err:
        out[name] = np.array(str(err))
np.savez(sys.argv[1], **out)
"""


def _run(tmp_path, limit_mb):
    script = tmp_path / "low_memory.py"
    script.write_text(SCRIPT)
    out = tmp_path / "out_{}.npz".format(limit_mb)
    env = dict(os.environ)
    env.pop("TIGRE_SYCL_MEMORY_LIMIT_MB", None)
    if limit_mb is not None:
        env["TIGRE_SYCL_MEMORY_LIMIT_MB"] = str(limit_mb)
    subprocess.run([sys.executable, str(script), str(out)], env=env, check=True, timeout=600)
    return dict(np.load(out))


@pytest.fixture(scope="module")
def unlimited(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("unlimited"), None)


def _assert_close(name, result, ref):
    assert result.dtype == np.float32, (name, str(result))
    assert np.isfinite(result).all(), name
    assert result.shape == ref.shape, name
    diff = np.abs(result - ref).max() / np.abs(ref).max()
    tolerance = {"AxSiddon": 1e-4}.get(name, 1e-5)
    assert diff < tolerance, (name, diff)


def test_split_results_match_the_unsplit_ones(tmp_path, unlimited):
    # 24 MB: Ax (image and projection buffers) splits.
    split = _run(tmp_path, 24)
    for name, ref in unlimited.items():
        _assert_close(name, split[name], ref)


def test_many_image_pieces(tmp_path, unlimited):
    # 2 MB: Ax splits the 64 slices into 11 pieces of 6 (12 pieces of 6 would leave the
    # last one with a negative number of slices).
    result = _run(tmp_path, 2)
    for name in ("Ax", "AxSiddon"):
        _assert_close(name, result[name], unlimited[name])


# A 16 MB image (64 slices, 256 kB each) and the small detector of SCRIPT.
# PROJ_PER_BLOCK is 9 in the cone projectors; one block is then ~0.14 MB.
_MEM_IMAGE = 64 * 256 * 256 * 4
_MEM_SLICE = 256 * 256 * 4
_MEM_PROJ_BLOCK = 9 * 61 * 67 * 4


def _expected_splits(mem_budget, max_alloc):
    # Same rule as tigreSyclImageSplits, then the no-empty-piece rounding the
    # projectors apply before allocating.
    n_voxel_z = 64
    if (
        _MEM_IMAGE + 2 * _MEM_PROJ_BLOCK < mem_budget
        and _MEM_IMAGE <= max_alloc
        and _MEM_PROJ_BLOCK <= max_alloc
    ):
        splits = 1
    else:
        if _MEM_SLICE > max_alloc or _MEM_PROJ_BLOCK > max_alloc:
            raise AssertionError("cannot even allocate one slice")
        if mem_budget <= 4 * _MEM_PROJ_BLOCK + _MEM_SLICE:
            raise AssertionError("budget too small for one slice")
        max_piece = min(mem_budget - 4 * _MEM_PROJ_BLOCK, max_alloc)
        splits = (_MEM_IMAGE + max_piece - 1) // max_piece
        splits = min(max(splits, 1), n_voxel_z)
    splits = (
        n_voxel_z
        + (n_voxel_z + splits - 1) // splits
        - 1
    ) // ((n_voxel_z + splits - 1) // splits)
    return splits


def _run_logged(tmp_path, name, limit_mb=None, max_alloc_mb=None):
    script = tmp_path / "low_memory.py"
    script.write_text(SCRIPT)
    out = tmp_path / "out_{}.npz".format(name)
    env = dict(os.environ)
    env.pop("TIGRE_SYCL_MEMORY_LIMIT_MB", None)
    env.pop("TIGRE_SYCL_MAX_ALLOC_MB", None)
    env["TIGRE_SYCL_SPLIT_LOG"] = "1"
    if limit_mb is not None:
        env["TIGRE_SYCL_MEMORY_LIMIT_MB"] = str(limit_mb)
    if max_alloc_mb is not None:
        env["TIGRE_SYCL_MAX_ALLOC_MB"] = str(max_alloc_mb)
    result = subprocess.run(
        [sys.executable, str(script), str(out)],
        env=env,
        check=True,
        timeout=600,
        capture_output=True,
        text=True,
    )
    logged = [
        int(line.split("=", 1)[1])
        for line in result.stderr.splitlines()
        if line.startswith("TIGRE_SYCL_SPLITS=")
    ]
    return dict(np.load(out)), logged


def test_max_alloc_cap_splits_without_over_splitting_the_budget(tmp_path, unlimited):
    # 8 MB is enough for half the 16 MB image as a single allocation, so two
    # pieces suffice. Treating 8 MB as the total budget (MEMORY_LIMIT_MB)
    # leaves less after the projection buffers and needs a third piece.
    split, logged = _run_logged(tmp_path, "maxalloc8", max_alloc_mb=8)
    expected = _expected_splits(mem_budget=(1 << 62), max_alloc=8 << 20)
    assert expected == 2, expected
    assert logged == [expected, expected], logged  # interpolated and Siddon
    # Interpolated matches the unsplit result. Two equal pieces put the split
    # exactly on the mid-plane of this centred cone geometry, so Siddon drops
    # the centre detector row in both pieces (the same mid-plane miss as in
    # the CUDA split path, #859). Check Siddon against a 3-piece budget cap.
    _assert_close("Ax", split["Ax"], unlimited["Ax"])

    budget, budget_logged = _run_logged(tmp_path, "limit8", limit_mb=8)
    for name, ref in unlimited.items():
        _assert_close(name, budget[name], ref)
    expected_budget = _expected_splits(mem_budget=8 << 20, max_alloc=(1 << 62))
    assert expected_budget > expected, (expected_budget, expected)
    assert budget_logged == [expected_budget, expected_budget], budget_logged


def test_large_max_alloc_does_not_split_when_the_budget_fits(tmp_path, unlimited):
    # A generous per-allocation cap must not force a split of a 16 MB image
    # that already fits in the total budget.
    result, logged = _run_logged(tmp_path, "maxalloc64", max_alloc_mb=64)
    for name, ref in unlimited.items():
        _assert_close(name, result[name], ref)
    assert logged == [1, 1], logged
