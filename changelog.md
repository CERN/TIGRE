# Changelog

All notable changes to TIGRE are listed here, newest first. Add a line under
**Unreleased** with every pull request that changes behaviour, and move them under
a version heading when releasing. Many contributors have made this possible:
thank you.

## Unreleased

### Added
- `from_vector_to_TIGRE` and `from_TIGRE_to_vector` (Python: `tigre.from_vector_to_TIGRE`,
  `tigre.from_TIGRE_to_vector`; MATLAB: `Utilities/`): convert between a cone-beam geometry and per-view
  vectors (source, detector centre, u, v) in the layout of ASTRA's `cone_vec`, with the first pixel as an
  option instead of the centre. Any rigid placement of source and detector relative to the image converts
  both ways, including static source arrays with no common rotation centre. Vector geometries TIGRE cannot
  describe (u and v not perpendicular, pixel size changing between views) are refused, not approximated.
  No new geometry fields and no kernel changes.

### Fixed
- Python geometry helpers in `tigre.utilities.common_geometry` now build the scans they describe; they are
  rebuilt on `from_vector_to_TIGRE`, with the same arguments:
  - `staticDetLinearSourceGeo`: the "fixed" detector moved by several mm between source positions, and the
    function failed on the default geometry (its integer `rotDetector`).
  - `ArbitrarySourceDetMoveGeo`: the "fixed" detector moved by several mm, sources spread around the object
    gave geometries hundreds of mm off (or NaN), and bad input raised `TypeError` instead of a message.
    `d_rot` is still added to `rotDetector`, in its (roll, pitch, yaw) order.
  - `ArbitrarySourceDetectorFixedObject`: detector offsets were scaled by `2*DSO/DSD` (right only at
    magnification 2), an offset detector also misplaced the source, and a source on the -x axis raised
    "Rotation matrix must be right handed". It no longer prints.
  - `rotation_from_vecs` returned a reflection (`-I`) for opposite vectors; it now returns a rotation.

## 3.3.0 (2026-11-21)

### Added
- `TIGRE_CUDA_ARCH` (Python build and MATLAB `Compile`): compile only for the given GPU architectures (e.g.
  `89`, `75,89`, or `native`), several times faster, for development and testing. Default unchanged: every
  architecture nvcc supports.
- Python and MATLAB builds: the GPU backend (`cuda`, `hip`, `sycl` or `metal`) is chosen with the
  `TIGRE_BACKEND` environment variable, or in MATLAB with `Compile -backend <name>`. When neither
  is given, the build uses the first backend whose compiler is installed, CUDA first, and prints
  which one it chose. Only `cuda` can be built so far; the
  others stop with a message pointing to the work that adds them (#747, #885, #876).
- Python: `tigre.Ax` and `tigre.Atb` take PyTorch CUDA tensors and CuPy arrays, and return
  the same kind of array on the same GPU, without copying through host memory. This uses the
  array's GPU only; if `gpuids` names other GPUs, the data goes through the host as before.
- PyTorch bindings: CUDA tensors stay on the GPU (see above). If TIGRE runs out of GPU memory,
  the bindings empty PyTorch's memory cache and try once more.
- Python demo `d25_Pytorch.py`: a section on passing PyTorch GPU tensors to `tigre.Ax`/`tigre.Atb`
  directly, with a timing against the round trip through the host; the demo now runs on GPU 0 and downloads
  MNIST to `./data/mnist` (it used GPU 3 and a path on one machine).
- Parallel beam: `Ax` (Siddon and interpolated) and `Atb` use all the GPUs in `gpuids`, as cone beam
  does: `Ax` shares the angles between the GPUs, `Atb` the image slices. Before, they ran on the first
  GPU only. Like cone beam, they now also work when the image does not fit on the GPU: they split it
  in pieces along z (before, they failed with a CUDA error) (#286, #295).

### Changed
- Backprojection (`Atb`, cone beam) is faster: the kernels step each voxel's position up its column instead of
  rebuilding it from the indices, the per-projection set-up is precomputed on the host, and `matched` computes its
  weight from numbers the kernel already has. Kernel time, RTX 4070: `matched` 0.73 of before, FDK weights 0.91;
  Quadro RTX 6000, measured before the host precomputation: 0.65 and 0.86. Results change by float rounding
  (~5e-5 of the maximum), and `matched` as described under Fixed.
- `Ax` and `Atb` with NumPy arrays no longer check `np.isreal` on the input, which the `float32` check already
  implies: it allocated two arrays the size of the input on every call (~0.1 s per call at 512³).
- `minTV`/`AwminTV` (the TV step of ASD_POCS, OS_ASD_POCS, AwASD_POCS, ...): the first and last
  z-slice of the image now use the same boundary as x and y (no neighbour beyond the edge). Before, the slices
  beyond the image's ends were read as voxels of value 0, so even a constant image was changed at its first and
  last slice. Results change there, and a little further in as the iterations go on (as #863 did for tvdenoise).
  `minPICCS` already treated z like x and y.
- `minPICCS` (PICCS) runs on the GPU code it now shares with `minTV`: it uses all the GPUs in `gpuids` (before,
  whichever GPU was current), works on images larger than the GPU memory, and sums the gradient norm in double.
  `PICCS.cu` keeps only the PICCS maths (the step direction); the multi-GPU and splitting code is
  `tv_gradient_descent` (`tv_gradient_descent.hpp`).
- `AwminTV` runs on the same shared code (`GD_AwTV.cu` keeps only the adaptive-weighted TV step
  direction): one copy of the multi-GPU and splitting code instead of two. Results are unchanged. Its error and
  warning identifiers are now `minTV`'s (`GD_TV:...`).
- `CTnoise.add` (Python) and `addCTnoise` (MATLAB) add a new noise realisation on every call. Before, the
  GPU random number generator was always seeded with 1234, so every call added exactly the same noise. The
  new seed is drawn from NumPy's global random state (Python) or MATLAB's random stream, so
  `np.random.seed(...)` / `rng(...)` make a script reproducible; `seed=` (Python) or `'Seed'` (MATLAB) fixes it.
- Python tests: rare edge cases and stress tests (little free GPU memory, thousands of repeated
  calls, images split between pieces) are marked `extended` and not run by a plain `pytest`
  (on an RTX 4070 the suite goes from ~6 to ~2 minutes). Run them with `pytest -m extended`;
  see `contributing.md`.
- `tvdenoise`/`im3DDenoise` (and so FISTA, SART-TV and the other TV algorithms that use it): the
  first and last z-slice of the image now use the same boundary as x and y. Before, a mirrored
  extra slice was added beyond each end and updated with the rest, so results at the first and
  last slice change (by up to ~10% of the maximum in strong denoising; under 1% one slice in;
  unchanged further in). This also makes split, unsplit and multi-GPU results identical (#863).
- `minTV`/`AwminTV` (used by ASD_POCS, OS_ASD_POCS, AwASD_POCS, ...): each GPU's part of the
  gradient norm is now summed and returned in double precision, and the parts are added in double,
  so they give the same result on 2 or more GPUs as on 1 (before, up to ~5e-4 of the maximum apart
  in one call). Single-GPU results change by rounding; the POCS solvers amplify such changes (up to
  ~5e-3 of the maximum after 10 OS_AwASD_POCS iterations) (#856; approach by Jeff Daily).
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
- `Atb` with `matched` weights used the wrong weight for rotated detectors (`rotDetector`) and for rotation axes
  that are not z (3-angle `angles`): it rebuilt the detector point as if the detector were not rotated and the
  source at z=0, and half a pixel off. The weight was off by 1.5-1.8% for a 2-degree roll or pitch, ~6.5% for
  5 degrees on all three or a 10-degree axis tilt, and 2e-4 for an unrotated detector. It is now exact for any
  geometry, and `matched` is as close to the adjoint of `Ax` for rotated detectors as for unrotated ones.
  `matched` results change: ~1e-4 for unrotated detectors, up to several % for rotated ones.
- `minPICCS` returned NaN when the gradient was zero (e.g. a constant image, or one PICCS had already
  flattened): it divided by the gradient norm without the small constant minTV adds. It now uses minTV's
  normalisation.
- MATLAB `minPICCS` with 5 inputs (no `gpuids`) left alpha, the iterations and the ratio unset, so it did
  nothing or something arbitrary; the prior's size and type were not checked (a double prior from `PICCS` was read
  as garbage). Now all input forms are read, the prior must match the image, and `minimizeTV_prior` passes singles.
- `AwminTV` (and `minTV`, through the shared code) on 3 or more GPUs, or split for memory, when the last piece was
  much shorter than the others: the gradient kernel's grid was sized from the last piece, so the other pieces
  could miss their last slices.
- Cone-beam `Ax`/`Atb` with a host array that is already page-locked (e.g. the NumPy view of a pinned PyTorch
  tensor) failed with "Error pinning memory", and unregistered memory that was not theirs. So did parallel-beam
  `Atb` and the TV functions with an array allocated page-locked (CUDA says "invalid argument" when asked to
  register it). All of them now use such an array as it is and leave it alone.
- `CTnoise.add` / `addCTnoise` ran on the current GPU whatever `gpuids` said, and leaked GPU memory if a CUDA
  error was raised (Python). They now run on the first of `gpuids` and free everything on error.
- `minTV`, `AwminTV` and `tvdenoise` failed on several GPUs (or when split for memory) when the image
  had so few z slices that the last pieces came out empty, e.g. 5 slices on 4 GPUs, or tvdenoise of
  one slice on 2 GPUs. The pieces are now chosen so that none is empty, using fewer GPUs if needed.
- `minTV`, `AwminTV` and `minPICCS` read past the end, or before the start, of the image, the output and
  PICCS's prior when a piece of the image (one per GPU, or more when it is split for memory) was thinner
  than the halo slices copied around its neighbour: e.g. 7 slices on 3 GPUs (pieces of 3, 3 and 1), or
  pieces of 1 slice. Depending on the memory there, this crashed, gave a copy error or did nothing visible;
  results were not affected. The halos now stop at the image's first and last slice, as in `tvdenoise`
  (#899; found by Einse57).
- `minTV`, `AwminTV`, `minPICCS` and `tvdenoise` freed some of their host memory only when they failed:
  every call lost a little, including a page-locked (pinned) allocation in `minTV`, `AwminTV` and `minPICCS`,
  which the POCS solvers call at every iteration.
- `tvdenoise` (Python; `im3DDenoise(...,'TV',...)` in MATLAB) on 2 or more GPUs without a direct link
  between them, when the image fits only shared between the GPUs: the host buffer for the halo slices
  had room for one slice per halo, but the halos are up to 60 slices wide, so they overwrote each other
  (wrong results) and the last one wrote past the end of the buffer (#891; found by Einse57 in #885).
- Python: `Ax` and `Atb` changed the `GpuIds` they were given when there were fewer angles
  (`Ax`) or image slices (`Atb`) than GPUs, keeping only as many GPUs as needed. Algorithms
  that project one angle at a time (SART, ASD_POCS, AwASD_POCS) then ran on one GPU after the
  first angle, for the projections and for the TV step. They now trim a copy (#855).
- `tvdenoise` (Python; `im3DDenoise(...,'TV',...)` in MATLAB, used by FISTA and SART-TV), when the
  image did not fit on the GPU and was split in pieces: for the last piece it zeroed a dual variable
  instead of the input; it did not wait for each piece's results to reach host memory before the next
  piece read from it; and a piece read its halo slices from host memory that the previous piece had
  already advanced. The result differed from the unsplit one across the whole volume (e.g. 7.6e-2 of
  the maximum with 100 MB free, up to ~2e-1 next to the joins). It now gives the unsplit result exactly
  (#806, #863).
- Cone-beam Siddon `Ax` got some rays wrong (by up to ~20%): for rays entering the image
  exactly on a voxel boundary, the first voxel and the first boundary crossings disagreed.
  Without splitting this was rare (e.g. 1 in 50,000 rays, 256³, off by up to 14%); when the
  image was split in pieces because it did not fit on the GPU, it was common. The traversal
  now derives both from the same voxel (#859).
- Siddon `Ax` (cone and parallel beam) is faster and more accurate. Cone beam, kernel time:
  1.4-1.9x faster on Windows (RTX 4070), 5.5-6.5x on Linux (Quadro RTX 6000), where the voxel
  indices were 64-bit integers. Faster: the loop fetches the voxel value once per step for all
  threads, instead of inside each branch, with 32-bit indices and no integer-to-float
  conversions. More accurate: it added up the step between voxel boundaries in single
  precision, which drifted by up to ~2e-3 (relative) along rays through a 512³ image. It now
  computes each boundary from its plane, with an error of ~1e-5. Siddon projections change by
  up to that drift, towards the exact value. Parallel beam uses the same traversal (#859), and is
  1.7-2.0x faster on Windows and 4.4-4.8x on Linux.
  The technical note `docs/technical-notes/siddon/siddon_note.pdf` explains these changes.
- `tvdenoise` (Python; `im3DDenoise` in MATLAB, also used by FISTA and other TV algorithms)
  could loop forever when the image did not fit on the GPU and was split: for some amounts of
  free memory it computed a buffer of 0 slices. `minTV`/`AwminTV` stopped with "Assertion
  Failed. Logic behind splitting flawed!" when free memory was low. All three now split
  correctly whenever the smallest piece fits, and raise "out of memory" otherwise (#858).
- `minTV`, `AwminTV`, `tvdenoise` and `minPICCS` now release everything they hold when they
  stop with an error (e.g. out of memory), as `Ax`/`Atb` do since #857. Before, an error left
  the caller's arrays registered as pinned memory; on Linux the next call, with a new array at
  the same address, then returned all zeros without any error, or failed with "part or all of
  the requested memory range is already mapped". A failed GPU allocation in these functions now
  reports "out of memory" instead of "invalid argument".
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
- Large problems: several sizes were computed with 32-bit integers and wrapped around. Matched
  `Atb` failed or left part of the image unscaled above 2^31 voxels on one GPU (e.g. 1300³);
  parallel-beam `Atb` could not allocate an image above 2^31 voxels. On Windows, where `long`
  is 32-bit: Python `Ax` allocated too small a buffer once the projections reached 16 GiB
  (e.g. 2048² x 1024 angles) and wrote past its end; `minTV`, `AwminTV`, `tvdenoise` and
  `minPICCS` got their sizes and offsets wrong above 2^31 voxels. Above ~15,400² detector pixels,
  `Ax` copy sizes overflowed too. All now use 64-bit sizes.

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
