"""CPU tests for the PyTorch bindings.

tigre.Ax and tigre.Atb are replaced by a small dense matrix and its transpose,
so these tests check the autograd plumbing, not the CUDA projectors.
"""

import numpy as np
import pytest

torch = pytest.importorskip("torch")

import tigre  # noqa: E402
from tigre.utilities import pytorch_bindings  # noqa: E402


@pytest.fixture
def operator(monkeypatch):
    geo = tigre.geometry_default(nVoxel=np.array([2, 3, 4]))
    geo.nDetector = np.array([2, 3])
    geo.dDetector = geo.sDetector / geo.nDetector
    angles = np.linspace(0, np.pi, 5, endpoint=False)
    n_img = int(np.prod(geo.nVoxel))
    n_proj = len(angles) * int(np.prod(geo.nDetector))
    M = np.random.default_rng(0).standard_normal((n_proj, n_img)).astype(np.float32)
    calls = []

    def Ax(img, geo, angles, *args, **kwargs):
        return (M @ img.ravel()).reshape(len(angles), *geo.nDetector).astype(np.float32)

    def Atb(proj, geo, angles, backprojection_type="FDK", **kwargs):
        calls.append(backprojection_type)
        return (M.T @ proj.ravel()).reshape(*geo.nVoxel).astype(np.float32)

    monkeypatch.setattr(tigre, "Ax", Ax)
    monkeypatch.setattr(tigre, "Atb", Atb)
    return geo, angles, M, calls


def test_shapes_with_batch_dimensions(operator):
    geo, angles, M, calls = operator
    A, At = pytorch_bindings.create_pytorch_operator(geo, angles, gpuids=["0"])
    x = torch.rand(2, 3, *geo.nVoxel.tolist())
    y = A(x)
    assert tuple(y.shape) == (2, 3, len(angles), *geo.nDetector.tolist())
    assert tuple(At(y).shape) == tuple(x.shape)


def test_gradient_of_A_uses_matched_backprojection(operator):
    geo, angles, M, calls = operator
    A, _ = pytorch_bindings.create_pytorch_operator(geo, angles, gpuids=["0"])
    x = torch.rand(1, *geo.nVoxel.tolist(), requires_grad=True)
    w = torch.rand(1, len(angles), *geo.nDetector.tolist())
    (A(x) * w).sum().backward()
    assert calls == ["matched"]
    expected = (M.T @ w.numpy().ravel()).reshape(x.shape)
    np.testing.assert_allclose(x.grad.numpy(), expected, rtol=1e-4, atol=1e-4)


def test_At_and_its_gradient(operator):
    geo, angles, M, calls = operator
    _, At = pytorch_bindings.create_pytorch_operator(geo, angles, gpuids=["0"])
    y = torch.rand(1, len(angles), *geo.nDetector.tolist(), requires_grad=True)
    w = torch.rand(1, *geo.nVoxel.tolist())
    (At(y) * w).sum().backward()
    assert calls == ["matched"]
    expected = (M @ w.numpy().ravel()).reshape(y.shape)
    np.testing.assert_allclose(y.grad.numpy(), expected, rtol=1e-4, atol=1e-4)


def test_fdk_weights_can_be_requested(operator):
    geo, angles, M, calls = operator
    A, At = pytorch_bindings.create_pytorch_operator(
        geo, angles, gpuids=["0"], backprojection_type="FDK"
    )
    x = torch.rand(1, *geo.nVoxel.tolist(), requires_grad=True)
    At(A(x)).sum().backward()
    assert calls == ["FDK", "FDK"]
