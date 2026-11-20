"""Ax and Atb with a host array that is already page-locked.

The operators page-lock (register) the caller's array to copy it faster. An
array that is already page-locked, e.g. the NumPy view of a pinned PyTorch
tensor, used to make the cone-beam operators fail with "Error pinning memory",
and at the end they unregistered memory that was not theirs. Now they use it as
it is and leave it alone.

Extended (needs PyTorch): pytest -m extended
"""
import numpy as np
import pytest

import tigre

pytestmark = [pytest.mark.gpu, pytest.mark.extended]
torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("PyTorch cannot use the GPU", allow_module_level=True)


def _pinned_copy(a):
    t = torch.empty(a.shape, dtype=torch.float32, pin_memory=True)
    view = t.numpy()
    view[...] = a
    return t, view


@pytest.mark.parametrize("mode", ["cone", "parallel"])
@pytest.mark.parametrize("op,kind", [("Atb", "FDK"), ("Atb", "matched"), ("Ax", "Siddon"), ("Ax", "interpolated")])
def test_pinned_input_gives_the_same_result(mode, op, kind):
    if mode == "parallel" and kind == "matched":
        pytest.skip("parallel beam has one backprojector")
    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array([32, 32, 32]))
    angles = np.linspace(0, 2 * np.pi, 20, endpoint=False)
    img = np.random.default_rng(0).random((32, 32, 32), dtype=np.float32)
    data = img if op == "Ax" else tigre.Ax(img, geo, angles)
    call = (lambda x: tigre.Ax(x, geo, angles, kind)) if op == "Ax" else (lambda x: tigre.Atb(x, geo, angles, kind))
    ref = call(data)
    t, view = _pinned_copy(data)
    for _ in range(2):  # the second call would fail if the first had unregistered the array
        np.testing.assert_array_equal(call(view), ref)
    assert t.is_pinned()
