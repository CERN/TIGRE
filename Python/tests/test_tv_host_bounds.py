"""minTV, AwminTV, minPICCS and tvdenoise read only the caller's arrays, never beyond them (#899).

These functions cut the image in z into pieces (one per GPU, and more if it does not fit), and copy each piece to
the GPU with up to buffer_length halo slices on each side. minTV, AwminTV and minPICCS copied buffer_length halo
slices on both sides of every piece but the first and the last, so next to a piece thinner than that the copy read
past the end of the image (e.g. 7 slices on 3 GPUs: pieces of 3, 3 and 1) or before its start (1-slice pieces).
That read the caller's image, the output and PICCS's prior out of bounds: a crash, a copy error or nothing visible,
depending on what lies there. The results were not affected. tvdenoise already limited its halos to the image.

Here the image (and the prior) sit between two pages that cannot be read, so any read outside them crashes the
process that runs the call (a subprocess, so that a crash fails the test).
- On 1 GPU, the GPU is filled (PyTorch) so that a 7 x 2048 x 2048 image is split into pieces of 2 slices (2, 2, 2,
  1) or of 1 slice; the result must still match the whole-GPU one (to ~1e-3: with several pieces per GPU these
  functions use a norm gathered from the previous pass, see test_tv_low_memory.py; tvdenoise to rounding).
- On 3 or more GPUs, few slices (4 up to 2 x GPUs + 1) are shared between all of them; the result must match 1 GPU.

Extended, Linux only (the protected pages use mmap/mprotect). Needs PyTorch for the 1-GPU case.
"""
import subprocess
import sys

import pytest

pytestmark = [pytest.mark.gpu, pytest.mark.extended]
if not sys.platform.startswith("linux"):
    pytest.skip("the protected pages use mmap/mprotect (Linux)", allow_module_level=True)

import tigre  # noqa: E402,F401  (first: CUDA's DLLs, and fail early if TIGRE is not built)
from tigre.utilities.gpu import GpuIds  # noqa: E402

try:
    from tigre.utilities.cuda_interface._minTV import minTV  # noqa: F401
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)
try:
    _N_GPUS = len(GpuIds().devices)
except Exception:
    _N_GPUS = 0
if _N_GPUS == 0:
    pytest.skip("no GPU", allow_module_level=True)

_COMMON = r"""
import ctypes, mmap, sys
import numpy as np
import tigre  # first
from tigre.utilities.gpu import GpuIds
from tigre.utilities.cuda_interface._minTV import minTV
from tigre.utilities.cuda_interface._AwminTV import AwminTV
from tigre.utilities.cuda_interface._tv_proximal import tvdenoise
from tigre.utilities.cuda_interface._minPICCS import minPICCS

_libc = ctypes.CDLL(None, use_errno=True)
_libc.mprotect.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
_keep = []


def guarded(values):
    # A copy of values between two pages that cannot be read or written.
    page = mmap.PAGESIZE
    nbytes = values.size * 4
    if nbytes % page:
        print("RESULT skip: the image is not a whole number of pages"); sys.exit(0)
    buf = mmap.mmap(-1, nbytes + 2 * page)
    addr = ctypes.addressof(ctypes.c_char.from_buffer(buf))
    for at in (addr, addr + page + nbytes):
        if _libc.mprotect(ctypes.c_void_p(at), page, 0) != 0:
            print("RESULT skip: mprotect failed"); sys.exit(0)
    _keep.append(buf)
    out = np.frombuffer(buf, dtype=np.float32, count=values.size, offset=page).reshape(values.shape)
    out[:] = values
    return out


def gpus(devices):
    ids = GpuIds()
    ids.devices = list(devices)
    return ids


def calls(ids):
    return {
        "minTV": lambda x, p: minTV(x, 0.01, 10, gpuids=ids),
        "AwminTV": lambda x, p: AwminTV(x, 0.01, 10, -0.005, gpuids=ids),
        "minPICCS": lambda x, p: minPICCS(x, p, 0.01, 10, 0.5, gpuids=ids),
        "tvdenoise": lambda x, p: tvdenoise(x, 10, 15.0, gpuids=ids),
    }
"""

# One GPU, split by filling the GPU. The image is 7 slices of 16 MiB. These functions use 95% of the free memory and
# need, per piece of s slices, n_images * (s + 4) slices (2 halo slices on each side): with n_images * (s + 4) <= the
# memory they use < n_images * (s + 5), they cut the image into pieces of s slices. n_images is 3 (minTV, AwminTV)
# or 4 (minPICCS: the prior too). tvdenoise splits as it needs.
_SPLIT = _COMMON + r"""
import torch

name, slices = sys.argv[1], int(sys.argv[2])
n_images = 4 if name == "minPICCS" else 3
slice_mib = 16
free = int((slices + 4.5) * n_images * slice_mib / 0.95 * 2**20)
rng = np.random.default_rng(0)
img = rng.random((7, 2048, 2048), dtype=np.float32)
prior = (0.5 * img).astype(np.float32)
call = calls(gpus([0]))[name]
reference = call(img.copy(), prior)


def fill(nbytes):
    tolerance = 8 * 2**20
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


x, p = guarded(img), guarded(prior)
ballast = fill(free)
if ballast is None:
    print("RESULT skip: could not leave {} MiB free".format(free >> 20)); sys.exit(0)
out = call(x, p)
del ballast
err = float(np.abs(out - reference).max() / np.abs(reference).max())
print("RESULT {:.3e}".format(err))
"""

# Several GPUs, few slices of 32 x 32 (4 KiB: one page).
_MULTI = _COMMON + r"""
name, nz, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(nz)
img = rng.random((nz, 32, 32), dtype=np.float32)
prior = (0.5 * img).astype(np.float32)
one = calls(gpus([0]))[name](img.copy(), prior)
every = calls(gpus(range(n)))[name](guarded(img), guarded(prior))
print("RESULT {:.3e}".format(float(np.abs(every - one).max() / np.abs(one).max())))
"""

NAMES = ["minTV", "AwminTV", "minPICCS", "tvdenoise"]


def _run(script, *args):
    try:
        r = subprocess.run([sys.executable, "-X", "faulthandler", "-c", script] + [str(a) for a in args],
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as e:
        pytest.fail("hung; output so far:\n{}".format(e.stdout))
    out = r.stdout + r.stderr
    line = next((ln for ln in r.stdout.splitlines() if ln.startswith("RESULT ")), None)
    assert r.returncode == 0 and line is not None, "crashed or failed (a read outside the arrays?):\n" + out
    line = line[len("RESULT "):]
    if line.startswith("skip"):
        pytest.skip(line)
    return float(line), out


@pytest.mark.parametrize("slices", [2, 1])  # pieces of 2, 2, 2, 1 slices: read past the end; of 1: before the start
@pytest.mark.parametrize("name", NAMES)
def test_split_on_one_gpu_reads_only_the_arrays(name, slices):
    pytest.importorskip("torch")
    import torch
    if not torch.cuda.is_available():
        pytest.skip("PyTorch cannot use the GPU")
    err, out = _run(_SPLIT, name, slices)
    assert err < (1e-6 if name == "tvdenoise" else 1e-3), out


@pytest.mark.skipif(_N_GPUS < 3, reason="needs 3 or more GPUs")
@pytest.mark.parametrize("nz", range(4, 2 * max(_N_GPUS, 3) + 2))
@pytest.mark.parametrize("name", NAMES)
def test_few_slices_on_all_gpus_read_only_the_arrays(name, nz):
    err, out = _run(_MULTI, name, nz, _N_GPUS)
    assert err < 1e-6, out
