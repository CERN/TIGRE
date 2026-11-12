# Projection operator contribution: contract and CGLS/SIRT

The operator owns geometry, GPU selection and masks. Missing measurements
are absent from data, residuals and normalization; masked NaN/Inf is removed
before arithmetic. Array layouts remain Python zyx/view-row-column and
MATLAB xyz/row-column-view. Subset indices are relative to the current
operator, zero based in Python and one based in MATLAB.

Conventional calls keep upstream defaults and implementations. Explicit
operator calls use selected-operator normalization and reject unsupported
geometry-specific options. Conventional Siddon/matched adapters are only
approximate adjoint pairs; exact transpose claims apply to the surface pair.

Python: `algs.cgls(data, operator=op, niter=20)` and `algs.sirt(...)`.
Use `TigreOperator(geo, angles)` or `SurfaceOperator(geo, surface, gpuids=0)`.
MATLAB: `CGLS_operator(data, op, 20)` / `SIRT_operator(...)`, or an explicit
`Operator` option on conventional names. Surface conveniences build an adapter
then call the same kernels. CGLS stops on zero normal residual/direction;
SIRT uses row and column sums, reciprocal zero on unobservable rays/voxels.
Operator defaults are zero initialization, relaxation 1, no nonnegativity;
MATLAB offers an optional relative residual tolerance (default 0).
Existing conventional defaults are not changed. FDK/FBP remain separate.
