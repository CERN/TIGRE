"""tvdenoise on 2 GPUs with a wide halo matches 1 GPU (#891).

When the image fits on the GPUs only shared between them (no further split),
but not 5 times on one, tvdenoise exchanges halos of buffer_length slices
between the GPUs, through host memory unless they have a direct link. The host
buffer had room for one slice per halo, so with buffer_length > 1 the halos
overwrote each other (wrong results) and the last one wrote past the end of
the buffer.

Here each GPU is left with 250 MB free for a 256^3 image (64 MB): 5 images do
not fit, half of them do, and the leftover gives buffer_length ~30. 100
iterations, so that halos are exchanged (every buffer_length iterations) and
used.

Extended (needs 2 GPUs and PyTorch, fills the GPUs): pytest -m extended
"""
import subprocess
import sys

import pytest

pytestmark = [pytest.mark.gpu, pytest.mark.extended]
torch = pytest.importorskip("torch")
if not torch.cuda.is_available() or torch.cuda.device_count() < 2:
    pytest.skip("needs 2 GPUs usable by PyTorch", allow_module_level=True)

FREE_MB = 250

_SCRIPT = r"""
import sys
import numpy as np
import torch
import tigre  # first: on Windows it puts CUDA's DLLs on the search path
from _tv_proximal import tvdenoise
from tigre.utilities.gpu import GpuIds

free_mb = int(sys.argv[1])
one = GpuIds(); one.devices = [0]
two = GpuIds(); two.devices = [0, 1]
img = np.random.default_rng(0).random((256, 256, 256), dtype=np.float32)
ref = tvdenoise(img.copy(), 100, 15.0, gpuids=one)
ballast = []
for d in (0, 1):
    with torch.cuda.device(d):
        torch.cuda.empty_cache()
        ballast.append(torch.empty(torch.cuda.mem_get_info(d)[0] - free_mb * 2**20, dtype=torch.uint8, device="cuda:%d" % d))
        torch.cuda.synchronize(d)
        if abs(torch.cuda.mem_get_info(d)[0] - free_mb * 2**20) > 16 * 2**20:
            print("RESULT could not leave that much free")
            sys.exit(0)
out = tvdenoise(img.copy(), 100, 15.0, gpuids=two)
print("RESULT {:.3e}".format(float(np.abs(out - ref).max() / np.abs(ref).max())))
"""


def test_two_gpus_wide_halo_match_one():
    try:
        r = subprocess.run([sys.executable, "-X", "faulthandler", "-c", _SCRIPT, str(FREE_MB)],
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as e:
        pytest.fail("hung; output so far:\n{}".format(e.stdout))
    out = r.stdout + r.stderr
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("RESULT ")), None)
    assert r.returncode == 0 and line is not None, out  # a heap overrun can also crash the process
    line = line[len("RESULT "):]
    if line == "could not leave that much free":
        pytest.skip(line)
    assert float(line) < 1e-6, out  # tvdenoise has no global sum: 2 GPUs give the 1-GPU result to rounding
