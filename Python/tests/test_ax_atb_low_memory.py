"""Ax and Atb with very little free GPU memory.

With 40 MB free for a 256^3 image (64 MB), cone-beam Ax must split the image in
pieces. Two problems showed up there (Linux, Quadro RTX 6000):

- The image texture's allocation (cudaMalloc3DArray) was not checked. When it
  failed, Siddon went on with an empty texture and the error it reported was
  "invalid resource handle" instead of "out of memory".
- The pieces were sized from the free memory as if a 3D array took exactly its
  raw size. It takes more (alignment, layout), so the allocation could fail
  where smaller pieces would have fitted.

Now every texture allocation reports "out of memory", and the forward
projectors allocate the texture first and split in more pieces if it does not
fit. Each case must either run and match the whole-GPU result, or raise a
clean "out of memory"; any other error fails.

The parallel-beam operators used to need the whole image on one GPU and failed
when it did not fit; now they split it too (#286, #295). Their cases use 251
slices, so that the last piece is shorter than the others.

Each case runs in its own process: on Windows the free memory CUDA reports
after a call on a nearly full GPU can stay low for a while, which would upset
the next case's fill.
"""
import subprocess
import sys

import pytest

pytestmark = pytest.mark.gpu
torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("PyTorch cannot use the GPU", allow_module_level=True)

FREE_MB = 40

_SCRIPT = r"""
import sys
import numpy as np
import torch
import tigre  # first: on Windows it puts CUDA's DLLs on the search path
from tigre.utilities.gpu import GpuIds

op, kind, free_mb, mode = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
ids = GpuIds()
ids.devices = [0]
nvoxel = [256, 256, 256] if mode == "cone" else [251, 256, 256]
geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array(nvoxel))
angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
img = np.random.default_rng(0).random(tuple(nvoxel), dtype=np.float32)
proj = tigre.Ax(img, geo, angles, "Siddon", gpuids=ids)
if op == "Ax":
    call = lambda: tigre.Ax(img, geo, angles, kind, gpuids=ids)
else:
    call = lambda: tigre.Atb(proj, geo, angles, kind, gpuids=ids)
reference = call()

torch.cuda.synchronize()
torch.cuda.empty_cache()
size = torch.cuda.mem_get_info()[0] - free_mb * 2**20
ballast = torch.empty(size, dtype=torch.uint8, device="cuda")
torch.cuda.synchronize()
now = torch.cuda.mem_get_info()[0]
if abs(now - free_mb * 2**20) > 8 * 2**20:
    print("RESULT could not leave that much free")
    sys.exit(0)
try:
    result = call()
    err = float(np.abs(result - reference).max() / np.abs(reference).max())
    print("RESULT ran {:.3e}".format(err))
except RuntimeError as e:
    print("RESULT raised {}".format(" ".join(str(e).split())))
"""

# The interpolated projector's split result differs slightly near piece borders (texture interpolation
# weights are 8-bit on the GPU; rays that skim a border have many samples there): ~1e-3 at most with the
# many small pieces 40 MB allows. Siddon and the backprojectors split to rounding.
CASES = [("Ax", "Siddon", 1e-4, "cone"), ("Ax", "interpolated", 3e-3, "cone"), ("Atb", "FDK", 1e-4, "cone"),
         ("Atb", "matched", 1e-4, "cone"),
         ("Ax", "Siddon", 1e-4, "parallel"), ("Ax", "interpolated", 3e-3, "parallel"), ("Atb", "FDK", 1e-4, "parallel")]


@pytest.mark.parametrize("op, kind, tolerance, mode", CASES)
def test_little_free_memory(op, kind, tolerance, mode):
    try:
        r = subprocess.run([sys.executable, "-c", _SCRIPT, op, kind, str(FREE_MB), mode],
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as e:
        pytest.fail("hung; output so far:\n{}".format(e.stdout))
    out = r.stdout + r.stderr
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("RESULT ")), None)
    assert r.returncode == 0 and line is not None, out
    line = line[len("RESULT "):]
    if line == "could not leave that much free":
        pytest.skip(line)
    if line.startswith("raised"):
        # What still fits in 40 MB depends on the platform and driver (CUDA's own
        # overheads), so a clean "out of memory" is acceptable; any other error is not.
        assert "out of memory" in line, out
        return
    assert line.startswith("ran"), out
    err = float(line.split()[1])
    assert err < tolerance, out
