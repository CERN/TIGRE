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
