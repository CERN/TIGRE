"""minTV, AwminTV and tvdenoise with little free GPU memory (#858).

When the image does not fit on the GPU, these split it in pieces along z. The
calculation of the pieces used unsigned sizes that could go below zero, and
tvdenoise could end up with a buffer of 0 slices and loop forever. On master a
256^3 image hung tvdenoise for most free-memory values between ~25 and ~170 MB,
and minTV/AwminTV stopped with "Assertion Failed. Logic behind splitting
flawed!" below ~25 MB.

Now: with enough memory for the smallest piece, each splits the image and
finishes; below that, each raises "out of memory" and works again afterwards.

What a split result should be:
- minTV/AwminTV: the whole-GPU result up to ~1e-3. They sum the gradient norm
  piece by piece in single precision (#856), so split results differ slightly.
- tvdenoise: only that it finishes with a finite result. Split, it differs from
  the whole-GPU result everywhere, by up to 26% here, on master too: a
  separate, older problem in how the pieces are combined.

The calls run in a separate process with a time limit, so a hang fails the test
instead of stopping the test run. The GPU is filled with a PyTorch tensor, and
the fill is checked and corrected after allocating (on Windows, other programs
drawing on the GPU make the free memory reported beforehand unreliable). If it
cannot leave the requested amount free, or another program takes memory during
a call, the case is skipped rather than failed.
"""
import subprocess
import sys

import pytest

pytestmark = pytest.mark.gpu
torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("PyTorch cannot use the GPU", allow_module_level=True)

# Free memory left for TIGRE, in MB: 2 is below the smallest piece (~4 MB for
# 256^3); the rest made master's tvdenoise hang or minTV/AwminTV fail.
FREE_MB = [2, 16, 30, 60, 116]

_SCRIPT = r"""
import sys
import numpy as np
import torch
import tigre  # first: on Windows it puts CUDA's DLLs on the search path
from _minTV import minTV
from _AwminTV import AwminTV
from _tv_proximal import tvdenoise
from tigre.utilities.gpu import GpuIds

ids = GpuIds()
ids.devices = [0]
img = np.random.default_rng(0).random((256, 256, 256), dtype=np.float32)
calls = {
    "minTV": lambda: minTV(img.copy(), 0.01, 10, gpuids=ids),
    "AwminTV": lambda: AwminTV(img.copy(), 0.01, 10, -0.005, gpuids=ids),
    "tvdenoise": lambda: tvdenoise(img.copy(), 10, 15.0, gpuids=ids),
}
reference = {name: call() for name, call in calls.items()}


def fill(nbytes):
    # Leave about nbytes free. Measure what is free after allocating and correct
    # the size: on Windows, with other programs drawing on the GPU, the free
    # memory reported before allocating can be off by a lot.
    tolerance = min(8 * 2**20, max(2 * 2**20, nbytes // 8))
    torch.cuda.synchronize()
    torch.cuda.empty_cache()
    size = torch.cuda.mem_get_info()[0] - nbytes
    for _ in range(4):
        if size <= 0:
            return None
        ballast = torch.empty(size, dtype=torch.uint8, device="cuda")
        torch.cuda.synchronize()
        now = torch.cuda.mem_get_info()[0]
        if abs(now - nbytes) <= tolerance:
            return ballast
        del ballast
        torch.cuda.empty_cache()
        size += now - nbytes
    return None


for free_mb in [int(a) for a in sys.argv[1:]]:
    for name, call in calls.items():
        ballast = fill(free_mb * 2**20)
        if ballast is None:
            print("{} {} MB | could not leave that much free".format(name, free_mb), flush=True)
            continue
        try:
            result = call()
            ref = reference[name]
            err = float(np.abs(result - ref).max() / np.abs(ref).max()) if np.isfinite(result).all() else float("inf")
            print("{} {} MB | ran, relative difference {:.2e}".format(name, free_mb, err), flush=True)
        except RuntimeError as e:
            torch.cuda.synchronize()
            if torch.cuda.mem_get_info()[0] < free_mb * 2**20 - 16 * 2**20:
                print("{} {} MB | another program took GPU memory".format(name, free_mb), flush=True)
            else:
                print("{} {} MB | raised, out of memory: {}".format(name, free_mb, "out of memory" in str(e)), flush=True)
        del ballast
        torch.cuda.empty_cache()
        np.testing.assert_array_equal(call(), reference[name])
        print("{} {} MB | works after".format(name, free_mb), flush=True)
"""


@pytest.mark.parametrize("free_mb", FREE_MB)
def test_tv_with_little_free_memory(free_mb):
    try:
        result = subprocess.run([sys.executable, "-X", "faulthandler", "-c", _SCRIPT, str(free_mb)],
                                capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as e:
        pytest.fail("hung (no result after 10 minutes); output so far:\n{}".format(e.stdout))
    out = result.stdout + result.stderr
    assert result.returncode == 0, out
    for name in ["minTV", "AwminTV", "tvdenoise"]:
        prefix = "{} {} MB | ".format(name, free_mb)
        line = next((ln for ln in result.stdout.splitlines() if ln.startswith(prefix)), None)
        assert line is not None, out
        if line.endswith(("could not leave that much free", "another program took GPU memory")):
            pytest.skip(line)  # other programs on the GPU changed its free memory under the test
        if free_mb == FREE_MB[0]:
            assert line == prefix + "raised, out of memory: True", out
        else:
            assert "ran" in line, out
            err = float(line.rsplit(" ", 1)[1])
            if name == "tvdenoise":
                assert err < float("inf"), out  # finite; see the module docstring
            else:
                assert err < 1e-3, out
        assert prefix + "works after" in result.stdout, out
