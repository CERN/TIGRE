"""minTV, AwminTV and tvdenoise on several GPUs with few z slices.

The TV functions cut the image in z into one piece per GPU (and per round, if
it does not fit), of ceil(slices / pieces) slices each. That can leave the last
pieces empty: 5 slices on 4 GPUs gives 2, 2, 1 and -1, and the -1 wraps around
(unsigned) to a huge size; an empty piece launches kernels with no blocks.
Either way the call failed. Now the pieces are chosen so that none is empty,
using fewer GPUs (or smaller pieces) when needed.

With 2 GPUs only tvdenoise of a single slice hit it (minTV and AwminTV use one
GPU below 4 slices); with 3 or more GPUs, many small sizes did.

Extended (needs 2 or more GPUs): pytest -m extended
"""
import numpy as np
import pytest

pytestmark = [pytest.mark.gpu, pytest.mark.extended]

import tigre  # noqa: E402,F401  (first: on Windows it puts CUDA's DLLs on the search path)
from tigre.utilities.gpu import GpuIds  # noqa: E402

try:
    from _minTV import minTV
    from _AwminTV import AwminTV
    from _tv_proximal import tvdenoise
    from _minPICCS import minPICCS
except ImportError:
    pytest.skip("TIGRE's CUDA extensions are not built", allow_module_level=True)

try:
    _ALL = list(GpuIds().devices)
except Exception:  # no GPU or no driver
    _ALL = []
if len(_ALL) < 2:
    pytest.skip("needs 2 or more GPUs", allow_module_level=True)


def _gpus(devices):
    ids = GpuIds()
    ids.devices = list(devices)
    return ids


CALLS = {
    "minTV": lambda img, ids: minTV(img, 0.01, 10, gpuids=ids),
    "AwminTV": lambda img, ids: AwminTV(img, 0.01, 10, -0.005, gpuids=ids),
    "tvdenoise": lambda img, ids: tvdenoise(img, 10, 15.0, gpuids=ids),
    "minPICCS": lambda img, ids: minPICCS(img, 0.5 * img, 0.01, 10, 0.5, gpuids=ids),
}


@pytest.mark.parametrize("nz", range(1, 2 * len(_ALL) + 2))
@pytest.mark.parametrize("name", sorted(CALLS))
def test_all_gpus_match_one_with_few_slices(name, nz):
    img = np.random.default_rng(nz).random((nz, 32, 40), dtype=np.float32)
    one = CALLS[name](img.copy(), _gpus(_ALL[:1]))
    every = CALLS[name](img.copy(), _gpus(_ALL))
    np.testing.assert_allclose(every, one, rtol=0, atol=1e-6 * float(np.abs(one).max()))
