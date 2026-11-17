"""tigre.Ax and tigre.Atb with GPU arrays (PyTorch CUDA tensors, CuPy arrays).

The data stays on the GPU, and the result must be exactly what the NumPy path
gives on the same GPU: the kernels are the same, only where the input comes
from and where the output goes changes.

Also checks the fallback (gpuids naming other GPUs: through the host, same
result), the checks on the input, that the producer's stream is waited for,
that nothing leaks, that a CUDA error raises and leaves TIGRE usable, and that
the PyTorch operators keep CUDA tensors on the GPU.
"""
import warnings

import numpy as np
import pytest

import tigre
from tigre.utilities import _gpu_arrays
from tigre.utilities.gpu import GpuIds

pytestmark = pytest.mark.gpu

torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("PyTorch cannot use the GPU", allow_module_level=True)
try:
    import cupy
except ImportError:
    cupy = None

FRAMEWORKS = ["torch"] + (["cupy"] if cupy is not None else [])


def to_gpu(a, framework, device=0):
    if framework == "torch":
        return torch.from_numpy(np.ascontiguousarray(a)).to("cuda:{}".format(device))
    with cupy.cuda.Device(device):
        return cupy.asarray(a)


def to_numpy(x):
    return x.cpu().numpy() if torch.is_tensor(x) else cupy.asnumpy(x)


def device_of(x):
    return x.device.index if torch.is_tensor(x) else x.device.id


def on_one_gpu(device=0):
    return GpuIds._from_devices([device])


def problem(mode, n=(24, 28, 32), nangles=40):
    geo = tigre.geometry(mode=mode, default=True, nVoxel=np.array(n))
    angles = np.linspace(0, 2 * np.pi, nangles, endpoint=False)
    rng = np.random.default_rng(0)
    img = rng.random(tuple(geo.nVoxel), dtype=np.float32)
    proj = tigre.Ax(img, geo, angles, gpuids=on_one_gpu())
    return geo, angles, img, proj


AX_CASES = [("cone", "Siddon"), ("cone", "interpolated"),
            ("parallel", "Siddon"), ("parallel", "interpolated")]
ATB_CASES = [("cone", "FDK"), ("cone", "matched"), ("parallel", "FDK")]


@pytest.mark.parametrize("framework", FRAMEWORKS)
@pytest.mark.parametrize("mode,projection_type", AX_CASES)
def test_ax_matches_numpy(framework, mode, projection_type):
    geo, angles, img, _ = problem(mode)
    expected = tigre.Ax(img, geo, angles, projection_type, gpuids=on_one_gpu())
    x = to_gpu(img, framework)
    result = tigre.Ax(x, geo, angles, projection_type)
    assert type(result) is type(x)
    assert device_of(result) == device_of(x)
    assert str(result.dtype).endswith("float32")
    np.testing.assert_array_equal(to_numpy(result), expected)
    np.testing.assert_array_equal(to_numpy(x), img)  # input untouched


@pytest.mark.parametrize("framework", FRAMEWORKS)
@pytest.mark.parametrize("mode,backprojection_type", ATB_CASES)
def test_atb_matches_numpy(framework, mode, backprojection_type):
    geo, angles, _, proj = problem(mode)
    expected = tigre.Atb(proj, geo, angles, backprojection_type, gpuids=on_one_gpu())
    y = to_gpu(proj, framework)
    result = tigre.Atb(y, geo, angles, backprojection_type)
    assert type(result) is type(y)
    assert device_of(result) == device_of(y)
    np.testing.assert_array_equal(to_numpy(result), expected)
    np.testing.assert_array_equal(to_numpy(y), proj)


def test_euler_angles():
    geo, _, img, _ = problem("cone")
    angles = np.random.default_rng(1).random((30, 3)) * 2 * np.pi
    expected = tigre.Ax(img, geo, angles, gpuids=on_one_gpu())
    result = tigre.Ax(to_gpu(img, "torch"), geo, angles)
    np.testing.assert_array_equal(result.cpu().numpy(), expected)
    np.testing.assert_array_equal(tigre.Atb(result, geo, angles, "matched").cpu().numpy(),
                                  tigre.Atb(expected, geo, angles, "matched", gpuids=on_one_gpu()))


def test_non_contiguous_input():
    geo, angles, img, _ = problem("cone")
    x = torch.from_numpy(np.ascontiguousarray(img.transpose(0, 2, 1))).cuda().transpose(1, 2)
    assert not x.is_contiguous()
    np.testing.assert_array_equal(tigre.Ax(x, geo, angles).cpu().numpy(),
                                  tigre.Ax(img, geo, angles, gpuids=on_one_gpu()))


def test_input_checks():
    geo, angles, img, proj = problem("cone")
    with pytest.raises(TypeError):
        tigre.Ax(to_gpu(img, "torch").double(), geo, angles)
    with pytest.raises(ValueError):
        tigre.Ax(to_gpu(img[:, :, :-1], "torch"), geo, angles)
    with pytest.raises(TypeError):
        tigre.Atb(to_gpu(proj, "torch").half(), geo, angles)
    with pytest.raises(ValueError):
        tigre.Atb(to_gpu(proj[:-1], "torch"), geo, angles)


def test_cpu_tensors_are_not_taken_for_gpu_arrays():
    assert not _gpu_arrays.is_gpu_array(torch.zeros(3))
    assert not _gpu_arrays.is_gpu_array(np.zeros(3))
    assert _gpu_arrays.is_gpu_array(torch.zeros(3, device="cuda"))


def test_explicit_gpuids_of_its_own_gpu_stay_on_the_gpu(monkeypatch):
    geo, angles, img, proj = problem("cone")

    def no_host(*args):
        raise AssertionError("went through the host")
    monkeypatch.setattr(_gpu_arrays, "_to_host", no_host)
    x = to_gpu(img, "torch")
    np.testing.assert_array_equal(tigre.Ax(x, geo, angles, gpuids=on_one_gpu()).cpu().numpy(), proj)
    tigre.Atb(tigre.Ax(x, geo, angles), geo, angles, gpuids=on_one_gpu())


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason="needs 2 GPUs")
def test_other_gpus_go_through_the_host_with_the_same_result():
    geo, angles, img, proj = problem("cone")
    x = to_gpu(img, "torch", device=1)
    both = GpuIds._from_devices([0, 1])
    result = tigre.Ax(x, geo, angles, gpuids=both)
    assert result.device.index == 1
    np.testing.assert_array_equal(result.cpu().numpy(), tigre.Ax(img, geo, angles, gpuids=both))
    # on its own GPU, the second one: GPU path on device 1
    result = tigre.Atb(to_gpu(proj, "torch", device=1), geo, angles)
    assert result.device.index == 1
    np.testing.assert_array_equal(
        result.cpu().numpy(), tigre.Atb(proj, geo, angles, gpuids=on_one_gpu(1)))


def test_waits_for_the_stream_that_wrote_the_input():
    """The input is written on a side stream, after a long kernel; TIGRE must
    not read it before that stream has finished."""
    geo, angles, img, _ = problem("cone")
    expected = tigre.Ax(img, geo, angles, gpuids=on_one_gpu())
    src = to_gpu(img, "torch")
    big = torch.rand(2048, 2048, device="cuda")
    side = torch.cuda.Stream()
    torch.cuda.synchronize()
    for _ in range(3):
        with torch.cuda.stream(side):
            x = torch.zeros_like(src)
            for _ in range(10):
                big = big @ big
                big /= big.max()
            x.copy_(src)
            result = tigre.Ax(x, geo, angles)
        side.synchronize()
        np.testing.assert_array_equal(result.cpu().numpy(), expected)


def test_repeated_calls_do_not_leak():
    geo, angles, img, proj = problem("cone", nangles=200)
    x, y = to_gpu(img, "torch"), to_gpu(proj, "torch")
    for _ in range(3):
        tigre.Atb(tigre.Ax(x, geo, angles), geo, angles)
    torch.cuda.synchronize()
    torch.cuda.empty_cache()
    before = torch.cuda.mem_get_info()[0]
    for _ in range(200):
        tigre.Ax(x, geo, angles)
        tigre.Atb(y, geo, angles, "matched")
    torch.cuda.synchronize()
    torch.cuda.empty_cache()
    drift_mb = (before - torch.cuda.mem_get_info()[0]) / 2**20
    assert drift_mb < 3, "{:.1f} MB of GPU memory lost over 400 calls".format(drift_mb)


def test_cuda_error_raises_and_tigre_keeps_working():
    geo, angles, img, proj = problem("cone")
    wide = tigre.geometry(mode="cone", default=True, nVoxel=np.array([1, 1, 131072]))
    with pytest.raises(RuntimeError):
        tigre.Ax(torch.ones(1, 1, 131072, device="cuda"), wide, angles)
    np.testing.assert_array_equal(tigre.Ax(to_gpu(img, "torch"), geo, angles).cpu().numpy(), proj)


def test_no_busy_gpu_warning_on_the_gpu_path():
    """Sharing the GPU with PyTorch is the point of the GPU path, so "the GPU is
    heavily used" is not reported there. The NumPy path still reports it."""
    free, total = torch.cuda.mem_get_info()
    grab = free - total // 2 + 256 * 2**20
    if grab <= 0 or grab > free - 512 * 2**20:
        pytest.skip("cannot leave the GPU between half full and nearly full")
    ballast = torch.empty(grab, dtype=torch.uint8, device="cuda")
    try:
        geo, angles, img, proj = problem("cone")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            tigre.Atb(to_gpu(proj, "torch"), geo, angles)
            tigre.Ax(to_gpu(img, "torch"), geo, angles)
        assert not [w for w in caught if "heavily used" in str(w.message)]
    finally:
        del ballast


# ---- the PyTorch operators
def test_pytorch_operators_keep_cuda_tensors_on_the_gpu(monkeypatch):
    from tigre.utilities.pytorch_bindings import create_pytorch_operator

    def no_host(*args):
        raise AssertionError("went through the host")
    monkeypatch.setattr(_gpu_arrays, "_to_host", no_host)

    geo, angles, img, proj = problem("cone")
    A, At = create_pytorch_operator(geo, angles)
    x = torch.from_numpy(np.stack([img, 2 * img])).cuda().requires_grad_()
    y = A(x)
    assert y.is_cuda and tuple(y.shape) == (2, len(angles), *geo.nDetector.tolist())
    np.testing.assert_array_equal(y[0].detach().cpu().numpy(), proj)

    w = torch.rand_like(y)
    (y * w).sum().backward()
    expected = np.stack([tigre.Atb(w[i].cpu().numpy(), geo, angles, "matched", gpuids=on_one_gpu())
                         for i in range(2)])
    np.testing.assert_array_equal(x.grad.cpu().numpy(), expected)

    single = At(torch.from_numpy(proj).cuda())
    assert single.is_cuda
    np.testing.assert_array_equal(single.cpu().numpy(),
                                  tigre.Atb(proj, geo, angles, "matched", gpuids=on_one_gpu()))
