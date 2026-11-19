# Explicit detector surfaces: model and implementation proposal

This proposes an opt-in explicit detector representation alongside conventional
TIGRE geometry. It is intended for modular flat panels, curved detector surfaces
and per-view pixel layouts that are awkward to describe with a single detector
plane. Existing geometry and projector calls remain the default.

This first patch adds the explicit model, a CUDA Siddon forward projector and
matching discrete adjoint, Python bindings, build integration and focused tests.
MATLAB bindings and opt-in shared solver support follow in dependent patches.
The existing native projector files and conventional Python calls are unchanged.

## Representation

`DetectorSurfaceGeometry` stores:

| Array | Shape | Meaning |
| --- | --- | --- |
| `source_positions_mm` | `(views, 3)` | Source positions in world `(x, y, z)`, mm |
| `pixel_centers_mm` | `(views, rows, columns, 3)` | Individual detector pixel centres in the same frame |
| `valid_mask` | `(views, rows, columns)` | Present measurements; false entries are absent rays |

A static `(rows, columns, 3)` pixel table or `(rows, columns)` mask may be
broadcast over the views. Stored geometry is finite contiguous float32 and the
mask is boolean. Geometry is separate from projection data. Volume-array order
remains TIGRE's existing convention; coordinate components here are `(x,y,z)`.

`DetectorModule` describes local pixel centres plus one rigid pose or one pose
per view. `from_modules` places modules in a common dense array and marks empty
layout cells invalid. Overlapping placements are rejected. Packing and unpacking
helpers associate each module's projection values with that layout. Gap values
in a packed array are placeholders; the surface projector/adjoint uses the mask
rather than interpret them as measured zero attenuation.

Subsets preserve the requested view order and own their arrays. Duplicate views
are allowed; empty, fractional, negative and out-of-range indices are rejected.
Pixel values may represent nonuniform spacing or genuinely curved surfaces;
no interpolation, planar fitting or gap completion is performed by this model.

## Example: two panels and an explicit gap

```python
import numpy as np
from tigre.utilities.detector_geometry import DetectorModule, DetectorSurfaceGeometry

pixels = np.zeros((2, 3, 3))
pixels[..., 1] = [-1., 0., 1.]
pixels[1, :, 2] = 1.
left_pose = np.eye(4)
left_pose[:3, 3] = [-100., -3., 0.]
right_pose = np.eye(4)
right_pose[:3, 3] = [-100., 3., 0.]
surface = DetectorSurfaceGeometry.from_modules(
    [[100., 0., 0.]],
    [DetectorModule(pixels, left_pose),
     DetectorModule(pixels, right_pose, column_offset=5)],
)
assert surface.valid_mask.shape == (1, 2, 8)
assert not surface.valid_mask[:, :, 3:5].any()
data = surface.pack_module_data([np.ones((1, 2, 3)), np.ones((1, 2, 3))])
```

## Projector API and numerical scope

```python
from types import SimpleNamespace
from tigre.utilities.Ax_surface import Ax_surface, Atb_surface

geo = SimpleNamespace(
    nVoxel=np.array([8, 16, 16]),       # z, y, x
    dVoxel=np.array([1., 1., 1.]),     # mm, z, y, x
    offOrigin=np.array([0., 0., 0.]),  # mm, z, y, x
)
volume = np.ones(tuple(geo.nVoxel), dtype=np.float32)
projections = Ax_surface(volume, geo, surface, gpuids=0)
backprojected = Atb_surface(projections, geo, surface, gpuids=0)
assert projections.shape == surface.valid_mask.shape
assert backprojected.shape == volume.shape
```

Each ray is the finite segment from its source to its detector pixel, with
Siddon voxel-intersection lengths in millimetres. This supports finite-source
cone/fan geometries with explicit per-view coordinates. Exact parallel beams
are outside this patch; a distant source is not an exact parallel-ray model.
No circular-trajectory assumption is used by these ray operators. That does
not establish analytical reconstruction validity for arbitrary trajectories.

Volume data and volume geometry fields use `(z, y, x)` order; explicit world
coordinates use `(x, y, z)`. Only `nVoxel`, `dVoxel` and `offOrigin` are read from
the conventional geometry object. The volume grid and origin must be constant
across views. Input volumes/projections are float32. Geometry is validated as
finite float32, and voxel spacing must be positive. Absent measurements produce
zero forward output and contribute nothing to the adjoint; the public adjoint
discards masked NaNs without modifying its input, while valid data must be finite.

Both directions use the same traversal weights, with the adjoint implemented
using float32 atomic accumulation. The transpose relation therefore holds up
to floating-point and summation-order errors; this is not FDK backprojection.
The traversal inherits TIGRE-style treatment of ray components smaller than
0.001 voxel as zero, so near-axis rays have that numerical approximation.

`gpuids` accepts a device index, indices, or `GpuIds`; the default selects all
visible devices. `n_slabs` partitions the volume along z, defaulting to one slab
per selected device. Each slab must fit its GPU. Forward slab contributions
are added, and adjoint slabs are concatenated. Detector tables are uploaded in
view chunks; this is a host NumPy API, not a GPU-resident CuPy pipeline.

The traversal follows [Siddon (1985)](https://doi.org/10.1118/1.595715),
also recorded in the Python utilities bibliography.

## Why explicit pixels and a shared projector contract?

Representable planar cases should continue to use conventional TIGRE geometry.
Masking conventional measurements is a valid way to remove detector samples.
Multiple planar modules may also be handled by composing per-panel conventional
projector calls, with an adjoint that combines their contributions.

The proposed explicit model provides a uniform description when pixel positions
are curved, irregular, or supplied per view. A mask changes which measurements
exist; it does not move their rays. We propose a surface projector for these
coordinates and a shared forward/adjoint/subset contract so iterative solvers
do not need separate permanent implementations for every geometry type.

| Route | Strength | Review question |
| --- | --- | --- |
| Conventional geometry plus masks | Reuses established kernels for representable planar cases | Which layouts fit the current geometry without changing the requested rays? |
| Per-panel conventional composition | Reuses kernels for multiple planar modules | What are the call, state, normalization and memory costs, and how are subsets composed? |
| Explicit pixels plus surface projector | One coordinate/mask model for planar and nonplanar rays | Does added native code justify itself through accuracy, maintenance and measured performance? |
| Virtual-detector rebinning | Enables a regular analytical reconstruction path | What interpolation and missing-data error does it introduce? |

These are design tradeoffs, not a measured claim that a new kernel always wins.
The review package should compare identical rays and masks on flat, tilted,
three-panel U and curved fixtures. Report coordinate agreement, analytical ray
integrals, adjointness, reconstruction error, memory and cold/warm runtime.

## Staged implementation

1. **CUDA/Python slice:** explicit model, rigid modules, layout masks, subsets,
   packing, native forward/adjoint pair, bindings and tests (this PR).
2. **MATLAB bindings:** reuse the same native pair; add MEX/build integration,
   independent array-layout checks and failure/recovery tests. Depends on PR 1.
3. **Python shared solver integration:** opt-in CGLS and SIRT with consistent
   masks/subsets and unchanged conventional call behavior. The intended operator
   contract is `forward(volume)`, `adjoint(data)` and `subset(view_indices)`, plus
   volume/data shapes and a validity mask. Its final public form is open to review.
   Depends on PR 1.
4. **MATLAB shared solver integration:** matching CGLS/SIRT contract and numerical
   checks, depending on the MATLAB bindings.
5. **Further families and utilities:** coherent solver groups and interchange
   after the initial model/backend contract has been reviewed.
6. **Analytical surface reconstruction:** a separate proposal with circular-orbit
   restrictions, rebinning and gap-completion policies made explicit. Do not infer
   that a geometry supported by ray projectors is valid for FDK.

The desktop application stays in a separate repository. Beam profiles,
collimation metadata, scan recipes, GUI exporters, animation, motion workflows,
CuPy acceleration and FDK implementations are outside this first patch.

## Cost and validation limits

Dense coordinates and a mask cost about `13 * views * rows * columns` bytes,
plus `12 * views` bytes for source coordinates. Constructor/broadcast and subset
copies can increase peak memory. Compact/factorized geometry and large-acquisition
limits should be considered before calling this a production-scale replacement.

Run the model tests without importing compiled TIGRE modules:

```shell
python -m pytest Python/tests/test_detector_surface_geometry.py -q
```

The model tests load its file directly because the existing package initializer
imports CUDA modules. After building and installing the package, run both suites:

```shell
python -m pytest Python/tests/test_detector_surface_geometry.py Python/tests/test_surface_projectors.py -q -rs
```

Native tests compare both directions with an independent float64 matrix formed
by intersecting each ray with every voxel box, plus axis-aligned uniform-box
integrals, flat agreement with conventional Siddon, masks, slab decomposition
and recovery after invalid-device errors. They do not demonstrate physical
detector response, measured-data reconstruction, multi-GPU scaling, absence of
all resource leaks, or an algorithmic speed advantage. The native suite skips
when `_AxSurface` or a CUDA GPU is unavailable; skipped tests are not native validation.
