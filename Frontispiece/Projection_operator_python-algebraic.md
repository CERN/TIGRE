# Projection operator contribution: algebraic subsets and TV

The operator owns geometry, GPU selection and masks. Missing measurements
are absent from data, residuals and normalization; masked NaN/Inf is removed
before arithmetic. Array layouts remain Python zyx/view-row-column and
MATLAB xyz/row-column-view. Subset indices are relative to the current
operator, zero based in Python and one based in MATLAB.

Conventional calls keep upstream defaults and implementations. Explicit
operator calls use selected-operator normalization and reject unsupported
geometry-specific options. Conventional Siddon/matched adapters are only
approximate adjoint pairs; exact transpose claims apply to the surface pair.

SART and OS-SART use subset-specific row/column normalization. Fast OS-SART
adds heuristic momentum; it does not guarantee monotone convergence.
TV regularization reuses existing volume kernels. Python upstream `ossart`
and `ossart_tv` remain available; `os_sart` / `os_sart_tv` are additive aliases.
Python fast OS-SART is operator-only; this slice does not introduce a new
conventional fast OS-SART implementation. SART and TV variants default to
nonnegativity on the operator path; operator OS-SART defaults to false.
