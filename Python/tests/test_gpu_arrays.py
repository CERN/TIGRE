"""tigre.Ax and tigre.Atb with GPU arrays (PyTorch CUDA tensors, CuPy arrays).

The data stays on the GPU, and the result must be exactly what the NumPy path
gives on the same GPU: the kernels are the same, only where the input comes
from and where the output goes changes.

The cases, by where the data is, how many GPUs, and whether it fits:

- GPU array, one GPU, fits: every Ax/Atb variant equals NumPy, bit for bit
  (test_ax_matches_numpy, test_atb_matches_numpy, Euler angles,
  non-contiguous input, input that requires grad).
- GPU array, GPU nearly full (other tensors, e.g. model weights), so TIGRE
  must split the image: same result (test_image_split_when_the_gpu_is_nearly_full).
- GPU array, gpuids naming other or several GPUs: through the host, same
  result (test_fallback_through_the_host; with 2 GPUs,
  test_other_gpus_go_through_the_host_with_the_same_result).
- PyTorch's cache holding the GPU: the operators still work
  (test_pytorch_operator_when_pytorch_cache_holds_the_gpu).
- The producer's stream is waited for; nothing leaks; a CUDA error raises and
  TIGRE keeps working (all 7 entry points); no "heavily used" warning.
- PyTorch operators: CUDA tensors stay on the GPU; CPU tensors go through NumPy.

NumPy input, the other half of the table, is covered by the existing tests
(test_projectors_gpu.py, test_cuda_errors.py, test_texture_object_leak.py):
that path is unchanged.
"""
import warnings

import numpy as np
import pytest

import tigre
from tigre.utilities import _gpu_arrays
from tigre.utilities.gpu import GpuIds

pytestmark = [
    pytest.mark.gpu,
    # The NumPy reference calls warn when PyTorch's cache holds over half of the
    # GPU, which is expected here; the test that checks the warning records it.
    pytest.mark.filterwarnings("ignore:.*heavily used:RuntimeWarning"),
]

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


def leave_free(nbytes):
    """Allocate a PyTorch tensor so that about nbytes of the GPU stay free."""
    torch.cuda.synchronize()
    torch.cuda.empty_cache()
    free, _ = torch.cuda.mem_get_info()
    if free <= nbytes + 64 * 2**20:
        pytest.skip("not enough free GPU memory for this test")
    try:
        return torch.empty(free - nbytes, dtype=torch.uint8, device="cuda")
    except RuntimeError:
        pytest.skip("could not fill the GPU")


def give_back_cache():
    torch.cuda.synchronize()
    torch.cuda.empty_cache()


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


def test_input_that_requires_grad():
    geo, angles, img, proj = problem("cone")
    x = to_gpu(img, "torch").requires_grad_()
    result = tigre.Ax(x, geo, angles)
    assert not result.requires_grad
    np.testing.assert_array_equal(result.cpu().numpy(), proj)


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


@pytest.mark.parametrize("framework", FRAMEWORKS)
def test_fallback_through_the_host(framework, monkeypatch):
    """What happens when gpuids names other GPUs, on one GPU: pretend the GPU
    path is not possible, so the data goes to the host and back."""
    geo, angles, img, proj = problem("cone")
    monkeypatch.setattr(_gpu_arrays, "_gpu_path_ids", lambda gpuids, device: None)
    result = tigre.Ax(to_gpu(img, framework), geo, angles, gpuids=on_one_gpu())
    assert type(result) is type(to_gpu(img, framework)) and device_of(result) == 0
    np.testing.assert_array_equal(to_numpy(result), proj)
    result = tigre.Atb(to_gpu(proj, framework), geo, angles, gpuids=on_one_gpu())
    np.testing.assert_array_equal(to_numpy(result), tigre.Atb(proj, geo, angles, gpuids=on_one_gpu()))


def test_which_gpuids_use_the_gpu_path():
    assert _gpu_arrays._gpu_path_ids(None, 0).devices == [0]
    assert _gpu_arrays._gpu_path_ids(on_one_gpu(0), 0).devices == [0]
    assert _gpu_arrays._gpu_path_ids([1], 1).devices == [1]
    assert _gpu_arrays._gpu_path_ids(on_one_gpu(1), 0) is None
    assert _gpu_arrays._gpu_path_ids(GpuIds._from_devices([0, 1]), 0) is None


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


TOO_WIDE = 131072  # wider than any GPU texture: a cheap, reliable CUDA error
ERROR_KINDS = [("Ax", "cone", "Siddon"), ("Ax", "cone", "interpolated"),
               ("Ax", "parallel", "Siddon"), ("Ax", "parallel", "interpolated"),
               ("Atb", "cone", "FDK"), ("Atb", "cone", "matched"), ("Atb", "parallel", "FDK")]


@pytest.mark.parametrize("op,mode,kind", ERROR_KINDS)
def test_cuda_error_raises_and_tigre_keeps_working(op, mode, kind):
    geo, angles, img, proj = problem("cone")
    bad_angles = angles[:4]
    if op == "Ax":
        wide = tigre.geometry(mode=mode, default=True, nVoxel=np.array([1, 1, TOO_WIDE]))
        with pytest.raises(RuntimeError):
            tigre.Ax(torch.ones(1, 1, TOO_WIDE, device="cuda"), wide, bad_angles, kind)
    else:
        wide = tigre.geometry(mode=mode, default=True, nVoxel=np.array([8, 8, 8]))
        wide.nDetector = np.array([1, TOO_WIDE])
        wide.dDetector = np.array([0.8, 0.8])
        wide.sDetector = wide.nDetector * wide.dDetector
        with pytest.raises(RuntimeError):
            tigre.Atb(torch.ones(4, 1, TOO_WIDE, device="cuda"), wide, bad_angles, kind)
    np.testing.assert_array_equal(tigre.Ax(to_gpu(img, "torch"), geo, angles).cpu().numpy(), proj)
    np.testing.assert_array_equal(tigre.Atb(to_gpu(proj, "torch"), geo, angles).cpu().numpy(),
                                  tigre.Atb(proj, geo, angles, gpuids=on_one_gpu()))


def test_no_busy_gpu_warning_on_the_gpu_path():
    """Sharing the GPU with PyTorch is the point of the GPU path, so "the GPU is
    heavily used" is not reported there. The NumPy path still reports it."""
    torch.cuda.empty_cache()
    free, total = torch.cuda.mem_get_info()
    ballast = None
    if free > total // 2 - 128 * 2**20:  # not busy yet: make it busy
        ballast = leave_free(total // 2 - 256 * 2**20)
    try:
        geo, angles, img, proj = problem("cone")
        x, y = to_gpu(img, "torch"), to_gpu(proj, "torch")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            tigre.Atb(y, geo, angles)
            tigre.Ax(x, geo, angles)
        assert not [w for w in caught if "heavily used" in str(w.message)]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            tigre.Atb(proj, geo, angles, gpuids=on_one_gpu())
        assert [w for w in caught if "heavily used" in str(w.message)], "NumPy path should still warn"
    finally:
        del ballast
        give_back_cache()


# ---- the GPU nearly full
SPLIT_CASES = [("Ax", "Siddon"), ("Ax", "interpolated"), ("Atb", "FDK"), ("Atb", "matched")]


@pytest.mark.parametrize("op,kind", SPLIT_CASES)
def test_image_split_when_the_gpu_is_nearly_full(op, kind):
    """Model weights or other tensors take most of the GPU: TIGRE cannot hold
    the whole image (64 MB) at once and works on it in pieces (only cone beam
    splits). About 42 MB stay free for TIGRE, well inside the range where it
    splits and each piece fits (10-70 MB for Ax, 20-80 MB for Atb)."""
    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([256, 256, 256]))
    angles = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    img = np.random.default_rng(0).random((256, 256, 256), dtype=np.float32)
    if op == "Ax":
        x = to_gpu(img, "torch")
        expected = tigre.Ax(img, geo, angles, kind, gpuids=on_one_gpu())
        out_bytes = expected.nbytes
    else:
        proj = tigre.Ax(img, geo, angles, gpuids=on_one_gpu())
        x = to_gpu(proj, "torch")
        expected = tigre.Atb(proj, geo, angles, kind, gpuids=on_one_gpu())
        out_bytes = img.nbytes
    ballast = leave_free(42 * 2**20 + out_bytes)  # the output is allocated first
    try:
        result = (tigre.Ax if op == "Ax" else tigre.Atb)(x, geo, angles, kind)
        torch.cuda.synchronize()
    finally:
        del ballast
        give_back_cache()
    # Splitting Ax changes the order of the sums: equal up to rounding.
    np.testing.assert_allclose(result.cpu().numpy(), expected, rtol=1e-5,
                               atol=1e-5 * float(np.abs(expected).max()))


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


def test_pytorch_operators_with_cpu_tensors():
    from tigre.utilities.pytorch_bindings import create_pytorch_operator

    geo, angles, img, proj = problem("cone")
    A, At = create_pytorch_operator(geo, angles)
    y = A(torch.from_numpy(img))
    assert not y.is_cuda
    np.testing.assert_array_equal(y.numpy(), tigre.Ax(img, geo, angles))
    np.testing.assert_array_equal(At(y).numpy(), tigre.Atb(proj, geo, angles, "matched"))


def test_pytorch_operator_when_pytorch_cache_holds_the_gpu(monkeypatch):
    """Training frees tensors into PyTorch's cache, which TIGRE cannot use.
    With all but ~3 MB of the GPU in that cache, TIGRE runs out of memory, the
    bindings give the cache back to CUDA and try again. (On Windows with the
    sysmem fallback on, CUDA may spill to system memory instead; the result
    must be right either way.)"""
    from tigre.utilities.pytorch_bindings import create_pytorch_operator

    geo = tigre.geometry(mode="cone", default=True, nVoxel=np.array([128, 128, 128]))
    geo.nDetector = np.array([512, 512])  # big projections: TIGRE needs >36 MB for them
    geo.dDetector = geo.sDetector / geo.nDetector
    angles = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    img = np.random.default_rng(0).random((128, 128, 128), dtype=np.float32)
    expected = tigre.Ax(img, geo, angles, gpuids=on_one_gpu())
    A, _ = create_pytorch_operator(geo, angles)
    x = to_gpu(img, "torch")

    emptied = []
    real_empty_cache = torch.cuda.empty_cache
    monkeypatch.setattr(torch.cuda, "empty_cache", lambda: (emptied.append(1), real_empty_cache()))
    ballast = leave_free(expected.nbytes + 3 * 2**20)
    emptied.clear()
    del ballast  # now in PyTorch's cache, not free for TIGRE
    try:
        result = A(x)
        torch.cuda.synchronize()
    finally:
        real_empty_cache()
    print("cache emptied {} time(s)".format(len(emptied)))
    np.testing.assert_array_equal(result.cpu().numpy(), expected)
