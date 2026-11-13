# Public API for external geometry-design clients

`tigre.design_api` 1.1 retains every 1.0 export and adds headless application
services. Recipes, explicit surfaces, YAML/NPZ/MAT interchange, physical
geometry checks and Matplotlib rendering remain in TIGRE. Desktop widgets,
event loops, dialogs and application state remain in the client package.
The extra `tigre[design]` installs PyYAML; importing TIGRE does not require YAML.

```python
from tigre import design_api as api

api.require_api(1, 1)
document = api.load_geometry_file("geometry.yaml")
# Saved counts/size/offsets must define the reconstruction volume. Convert
# moving-object geometry explicitly with fixed_object_geometry when needed.
api.validate_transmission_geometry(document.surface, document.geometry)
operator = api.create_operator(document.geometry, surface=document.surface,
                               gpuids=[0], n_slabs=2)
events = []
result = api.reconstruct(projections, operator=operator, algorithm="cgls",
                         niter=20, options={"return_residuals": True},
                         progress=events.append)
volume = result.volume
residual_history = result.diagnostics
```

`available_algorithms()` returns immutable descriptors for the registered
operator methods. Aliases `ossart` and `ossart_tv` remain accepted. Backend
registration is not proof of device availability or numerical convergence.
`available_devices()` explicitly queries zero-based CUDA device indices;
binding accepts the parent's device selections and surface slab options.
`backend_capabilities()` never queries CUDA and reports unsupported in-flight
cancellation and iteration progress explicitly. CPU-only imports do not imply
that the eager native TIGRE package can reconstruct without CUDA.

`create_operator` accepts exactly one of explicit `surface` or conventional
`angles`, owns its geometry/device settings and does not silently move,
resize, rebin or approximate the acquisition. Native and validation errors
remain actionable. Execution errors are wrapped in `BackendExecutionError`
with the original exception available through `__cause__`.

`reconstruct` preserves caller projection/initializer buffers and returns an
owned float32 volume, canonical algorithm name, requested iteration count,
solver-call duration and optional algorithm-specific diagnostics. Diagnostics
are a residual/error array, POCS history mapping, or None; no universal history
or actual completed-iteration count is invented. Supported options remain those
of the selected operator method. POCS requires an explicit absolute `maxl2err`.

Progress callbacks receive `validating`, `running` and `completed` lifecycle
events synchronously on the caller's thread. They do not run on the UI thread
automatically. Client applications marshal events through their own queues.
Callback exceptions propagate unchanged. Run reconstruction in a worker when
an interactive application needs a responsive event loop.

Cancellation is a callable predicate checked before validation, immediately
before execution and after execution. A running solver/native call completes
before cancellation can be observed. Its result is discarded if cancellation
is then requested. `ReconstructionCancelled` distinguishes this from a backend
failure. There is no per-iteration stop, partial result, checkpoint/resume or
promise of interrupting a CUDA kernel in this contribution.

Surface-FDK reconstruction is a separate contribution. To retain the existing
design API's remapping and geometry diagnostics without importing its deferred
native implementation, circular geometry checks live in a CPU helper. The
legacy `fdk_supported` field describes geometry restrictions only; the new
`fdk_available` field is false in this slice, as is the surface_fdk capability.
No experimental per-view reconstruction is enabled. Remapping remains the
existing explicit interpolation approximation, with physical panel gaps kept.
Rendering uses a local box helper rather than importing the fork's animation
rewrite. Existing native kernels, conventional projectors and solver
recurrences are unchanged.

This Python contribution depends on the prepared Python operator stack after
#814 is accepted. The campaign also waits for #815 acceptance and pending Linux
and real multi-GPU release checks. The separate MATLAB application still uses
its pinned named-helper contract. This Python API version does not certify a
MATLAB backend or replace its initialization/build requirements; a versioned
MATLAB facade and extraction of its design helpers remain a separate slice.
