# Projection operator contribution: ASD-POCS and PCSD

The operator owns geometry, GPU selection and masks. Missing measurements
are absent from data, residuals and normalization; masked NaN/Inf is removed
before arithmetic. Array layouts remain Python zyx/view-row-column and
MATLAB xyz/row-column-view. Subset indices are relative to the current
operator, zero based in Python and one based in MATLAB.

Conventional calls keep upstream defaults and implementations. Explicit
operator calls use selected-operator normalization and reject unsupported
geometry-specific options. Conventional Siddon/matched adapters are only
approximate adjoint pairs; exact transpose claims apply to the surface pair.

All eight ASD-POCS/PCSD variants require a finite nonnegative absolute masked
L2 tolerance (`maxl2err` / `MaxL2err`); no FDK-derived default is assumed.
PCSD uses residual-ratio TV scheduling and bounded iteration counts; legacy
mixed squared/unsquared gates are not imported. ASD and PCSD controls differ
and unsupported options fail. Weighted TV accepts scalar Delta only; vector
directional TV requires a separate native extension and is rejected here.
