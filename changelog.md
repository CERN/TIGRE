# Changelog

All notable changes to TIGRE are listed here, newest first. Add a line under
**Unreleased** with every pull request that changes behaviour, and move them under
a version heading when releasing. Many contributors have made this possible:
thank you.

## Unreleased

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
