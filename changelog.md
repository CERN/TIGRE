# Changelog

All notable changes to TIGRE are listed here, newest first. Add a line under
**Unreleased** with every pull request that changes behaviour, and move them under
a version heading when releasing. Many contributors have made this possible:
thank you.

## Unreleased

### Added
- Python: `tigre.Ax` and `tigre.Atb` take PyTorch CUDA tensors and CuPy arrays, and return
  the same kind of array on the same GPU, without copying through host memory. This uses the
  array's GPU only; if `gpuids` names other GPUs, the data goes through the host as before.
- PyTorch bindings: CUDA tensors stay on the GPU (see above). If TIGRE runs out of GPU memory,
  the bindings empty PyTorch's memory cache and try once more.

### Changed
- Python: a CUDA error now raises `RuntimeError` (with the CUDA message) instead of ending the
  whole Python process with `exit(1)`. Ax and Atb release the GPU memory, streams, textures and
  page-locked memory the failed call held, so TIGRE can be used again in the same process.
  In Python, TIGRE no longer calls `cudaDeviceReset()` after an error, which would destroy every
  other GPU allocation in the process (e.g. PyTorch's). MATLAB behaviour is unchanged (#853).
- "One (or more) of your GPUs is being heavily used by another program" (less than half of the
  GPU memory free) is now a warning instead of an error, in Python and MATLAB. In Python,
  TIGRE's warnings are now Python `RuntimeWarning`s, shown once, instead of text printed on every
  call.
- PyTorch bindings: `At` and the gradient of `A` now use the `matched` backprojector, the
  closest to the transpose of `Ax`; before, they used FDK weights, so gradients through
  TIGRE were further from the true ones. Pass `backprojection_type="FDK"` to
  `create_pytorch_operator` (or `A`/`At`) to get the old behaviour.
- Wheels are built for Python 3.11 to 3.14 (3.14 was missing from 3.2.0), on
  `manylinux_2_28`: they need glibc 2.28 or newer (e.g. RHEL/AlmaLinux 8, Ubuntu 18.10+);
  older Linux distributions should build from source.

### Fixed
- MATLAB `Compile`: compiles for the GPU architectures the installed `nvcc` supports, as it
  reports them (`nvcc --list-gpu-arch`), like the Python build. The fixed list asked CUDA
  12.2-12.7 for `compute_120`, which needs 12.8, and compilation failed.
- Cone-beam `Ax` and `Atb` with almost no free GPU memory (e.g. when PyTorch holds the GPU):
  `Atb` killed the process with an integer division by zero, and `Ax` tried the whole image
  and failed. Both now raise "out of memory" with how much memory they need, before
  allocating anything (Python: `RuntimeError`, which the PyTorch bindings retry after
  emptying PyTorch's cache; MATLAB: an error).
- `Ax` with `"interpolated"` (cone and parallel) and parallel-beam Siddon started projecting
  without waiting for the image to finish copying to the GPU, so the first block of 8-9
  projections could be computed from a partly copied image. Seen with GPU arrays; with NumPy
  input it could happen for cone-beam interpolated when the image was split (did not fit on
  the GPU), as the image is then page-locked and its copy fully asynchronous.
- `Atb` with `backprojection_type="matched"` (used by the Krylov algorithms, e.g. CGLS, LSQR)
  could return parts of the image as zeros when the image did not fit on the GPU in one piece:
  copying a piece out raced with clearing the GPU buffer for the next piece. Python and MATLAB.
- PyTorch bindings: `create_pytorch_operator(geo, angles)` without `gpuids`, and `gpuids` given
  as a list of device numbers (`["0"]`), failed: TIGRE needs a `GpuIds`. Both work now.
- Atb leaked texture objects on every call (all but two of every 32-angle block), until CUDA
  refused to create more, about a million later: "Texture object creation fail ... out of
  memory" after a few hundred iterations however much memory was free. Ax leaked one per image
  split (#617, #851).
- Python: `Atb`, `minTV`, `AwminTV`, `minPICCS`, `tvdenoise` and the noise generators no longer
  leak their small list of GPU ids on every call.
- `load_head_phantom` could not find `head.mat` after `pip install .` on Linux or after
  `pip install --user`: pip puts the file in `<prefix>/data`, which was not searched (#850).

## 3.2.0 (2026-11-11)

### Changed
- **Python 3.11 or newer is required**; Python 3.10 is no longer supported (#825).
- **TIGRE is no longer published to PyPI.** The `tigre` name there belongs to an unrelated
  package, so `pip install tigre` does *not* install TIGRE. Install from source or with
  `pip install git+https://github.com/CERN/TIGRE.git`; built wheels are attached to each
  GitHub Release (#819).
- FDK/FBP filtering in Python runs in cache-sized blocks on several threads (#790).
- FISTA chooses its step size automatically by default (power iteration) (#769).
- MATLAB SART, SIRT and PCSD are now thin wrappers over shared code (#786).
- Parker weighting detects near-360° scans more robustly (#743).
- FDK no longer applies Wang detector-offset weights on short scans (#778).

### Added
- Parallel-beam Siddon projector supports full Euler angles (alpha, theta, psi) (#792).
- MATLAB `plotSinogram` (#771).
- Continuous integration: the Python test suite runs on every push and pull request (#816).
- GPU test suite for the projectors and every reconstruction algorithm, run locally with
  `pytest -m gpu Python/tests` (#821).
- `CITATION.cff`, a rewritten contributing guide and a pull request template (#818, #817).

### Fixed
- `minTV`/`AwminTV` (used by ASD-POCS, PCSD, PICCS and the TV variants) could
  intermittently return an all-zero volume: a race between `cudaMemset` and stream
  work. The same pattern is fixed in the projectors and backprojectors (#823).
- Windows source builds failed with setuptools 84 or newer (CUDA files were sent to
  the C++ compiler) (#822).
- Multi-GPU TV halo exchange: synchronization and offsets (#804).
- Parallel-beam projections had stripe artifacts from wrong intersection bounds (#752).
- Siddon projector could hang on rays that miss the volume (#779).
- TV/PICCS minimizers destroyed the caller's CUDA context on success, and leaked
  GPU memory; IRN-TV-CGLS failed under NumPy 2 (#770, #800).
- PICCS: could hang computing W, and other Python PICCS fixes (#784, #794).
- POCS methods (ASD-POCS, PICCS, PCSD): `lambda_red` never reached the data step (#766).
- Krylov methods: float64 norm accumulation, working restarts, results no longer
  discarded (#774).
- Iterative algorithms: per-angle geometry indexing (MATLAB) and block size 1
  (Python) (#801, #742).
- Parker weights took the fan angle from the wrong axis (#767).
- `im3DNORM` returned `None` for supported norm orders (#805).
- Measure_Quality UQI now uses the N-1 variance, matching MATLAB (#796).
- FBP filter parameter was not propagated (#760).
- Geometry: single-angle input under NumPy 2, and `cast_to_single` turning integers
  into floats (#749, #750).
- `im3ddenoise` rejects 2D input with a clear error (#748).
- Varian loader rotation direction and XIM parsing (#762).
- `gpuids=None` now means all GPUs everywhere; helical demo offset axis; computeCOR
  port (#773).
- Importing a compiled module (e.g. `_minTV`) before `tigre` failed with a circular
  import; removed a deprecated `numpy.core` import (#824).
- Plotting: one-dimensional and ZYZ angle arrays, empty `savegif`, MATLAB quality
  histories (#812, #746, #813).
- Clearer error when the bundled head phantom is missing (#811).
- Demos: duplicate 2π projections, quality demo inputs, geometry demo header (#759, #810, #761).
- Build: compute capability table and detection, NumPy build pin for manylinux2014,
  documented GCC version (#732, #738, #798, #799).

### Removed
- Read the Docs configuration (the docs it pointed to no longer exist) and other stale
  files (#817).

## Earlier changes (since v2.2)

### Features

- Faster filtering for FDK in python
- Better timing info for python algorithms
- Bruker Skyscann Data Loader

### Fixed

- Fix L2 and quality measure arrays bug
- FDK detector offset weights missing a factor of 2
