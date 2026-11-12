# Projection operator contribution: ISTA/FISTA and MLEM

The operator owns geometry, GPU selection and masks. Missing measurements
are absent from data, residuals and normalization; masked NaN/Inf is removed
before arithmetic. Array layouts remain Python zyx/view-row-column and
MATLAB xyz/row-column-view. Subset indices are relative to the current
operator, zero based in Python and one based in MATLAB.

Conventional calls keep upstream defaults and implementations. Explicit
operator calls use selected-operator normalization and reject unsupported
geometry-specific options. Conventional Siddon/matched adapters are only
approximate adjoint pairs; exact transpose claims apply to the surface pair.

ISTA/FISTA use the selected forward/adjoint for every data step and power
estimate; TV is performed by existing native volume kernels. Operator MLEM
requires nonnegative valid data and uses masked sensitivity. Its denominator
floor is 1e-6 and multiplicative ray ratios are clipped to [0,10]. This is
an explicit stabilization policy, not unrestricted textbook EM or a claim
of likelihood monotonicity; conventional MLEM is unchanged.
