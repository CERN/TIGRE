# Projection operator contribution: Krylov and regularized Krylov

The operator owns geometry, GPU selection and masks. Missing measurements
are absent from data, residuals and normalization; masked NaN/Inf is removed
before arithmetic. Array layouts remain Python zyx/view-row-column and
MATLAB xyz/row-column-view. Subset indices are relative to the current
operator, zero based in Python and one based in MATLAB.

Conventional calls keep upstream defaults and implementations. Explicit
operator calls use selected-operator normalization and reject unsupported
geometry-specific options. Conventional Siddon/matched adapters are only
approximate adjoint pairs; exact transpose claims apply to the surface pair.

Reuse existing Krylov recurrences with per-instance projection bindings.
Include LSQR, LSMR, hybrid LSQR, IRN-TV-CGLS, flexible LSQR-TV and AB/BA-GMRES.
Exact-breakdown guards, scalar reorthogonalization, damped-LSMR stopping and
the flexible TV preconditioner are operator-only numerical policies.
The flexible preconditioner uses up to 50 mean-free CG steps, relative normal
residual target 1e-6. Flexible methods can stop early when their explicit
residual rises; this is not proof of convergence on arbitrary acquisitions.
FDK/multigrid initialization and FDK-based GMRES backprojection are rejected.
