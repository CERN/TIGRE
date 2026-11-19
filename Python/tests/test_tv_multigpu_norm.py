"""minTV and AwminTV give the same result on 2 GPUs as on 1 (#856).

Each step divides the TV gradient by its norm. Each GPU sums its part of the
norm and the host adds the parts. The parts used to come back in single
precision and be added in single precision, so the norm, and the whole
result, depended on how the image was shared between the GPUs: up to ~3.5e-3
of the image maximum for one call at the alpha the POCS solvers use, and
~2e-2 end to end for AwASD_POCS. Now the last step of each GPU's sum, and the
host's sum of the parts, are in double.

Needs 2 GPUs; skipped otherwise.
"""
import numpy as np
import pytest

pytestmark = pytest.mark.gpu

import tigre  # noqa: E402,F401  (first: on Windows it puts CUDA's DLLs on the search path)
from tigre.utilities.gpu import GpuIds  # noqa: E402

try:
    from _minTV import minTV
    from _AwminTV import AwminTV
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)


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


@pytest.mark.parametrize("n", [256, 250])  # 250: the GPUs' parts do not end on reduction-block boundaries
@pytest.mark.parametrize("alpha", [0.002, 0.6])
@pytest.mark.parametrize("name", ["minTV", "AwminTV"])
def test_two_gpus_match_one(name, alpha, n):
    img = np.random.default_rng(0).random((n, n, n), dtype=np.float32)
    if name == "minTV":
        call = lambda ids: minTV(img.copy(), alpha, 10, gpuids=ids)
    else:
        call = lambda ids: AwminTV(img.copy(), alpha, 10, -0.005, gpuids=ids)
    one = call(_gpus([0]))
    two = call(_gpus([0, 1]))
    err = float(np.abs(two - one).max() / np.abs(one).max())
    assert err < 1e-6, err
