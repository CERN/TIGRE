# Opt-in projection operators

`ProjectionOperator` supplies `forward`, `adjoint`, `subset`, `volume_shape`,
`projection_shape` and a boolean `valid_mask` (or None). Python arrays use zyx
volumes and view/row/column projections. Inputs are owned; missing NaN/Inf
samples are removed before arithmetic. Nested subsets preserve ordering.

Use `algs.cgls(data, operator=op, niter=20)` or `algs.sirt(...)` with a
`TigreOperator(geo, angles)` or `SurfaceOperator(geo, surface, gpuids=0)`.
`cgls_surface` and `sirt_surface` are thin conveniences over the same kernels.
The conventional adapter's Siddon/matched pair is approximate, not an exact
transpose; the surface adapter uses the discrete transpose from #814.

Conventional calls retain upstream implementations/defaults. Opt-in SIRT uses
selected-operator row/column sums, zero weights for unobservable samples/voxels,
relaxation 1 and nonnegativity disabled by default. CGLS accepts `init` and
`return_residuals`; SIRT accepts `init`, `lmbda`, `nonneg`. Geometry-specific
weights, initializers and unsupported options fail. FDK/FBP are not part of
this interface. GPU/device/slab settings belong to the adapter.

The surface adapter requires the explicit detector-surface projector extension.
MATLAB currently selects one device. Multi-GPU scalability is not qualified.
