# Versions and backwards compatibility

This page says what you can rely on when you upgrade TIGRE, and how things
are changed or removed. It applies to both the MATLAB and the Python
versions, which share one version number.

## Version numbers

TIGRE versions are `MAJOR.MINOR.PATCH` (for example 3.2.0):

- **Patch** releases fix bugs. They add no features and remove nothing.
- **Minor** releases add features and may deprecate things, but do not remove
  anything that was not already deprecated.
- **Major** releases may remove deprecated features and change conventions.

Every change that affects users is listed in the
[changelog](https://github.com/CERN/TIGRE/blob/master/changelog.md), with
removals under *Removed* and deprecations under *Deprecated*.

## What is covered

The **public interface** is what this documentation describes:

- the functions and options in the [API reference](../api/index.md), in
  both languages;
- the geometry: the `geo` fields, their units and meaning, and the axis and
  angle conventions described in [](geometry.md);
- the scanner loaders and what they return;
- the array layouts that TIGRE takes and returns.

Anything else is internal and may change in any release: modules and
functions whose names start with an underscore, the compiled CUDA modules
(`_Ax`, `_Atb`, ...), the MATLAB `Utilities/cuda_interface` files and the
demos.

## Geometry

Geometry descriptions are often written once and kept for years, sometimes
for a whole scanner, so they get the strongest guarantee. Within a major
version:

- **Field names, units and conventions do not change.** All lengths are in
  millimetres, angles in radians, and the coordinate system, the axis order
  of each language and the meaning of each field stay as documented.
- **New fields are optional.** When a field is added, its default
  reproduces the previous behaviour, so an existing geometry gives the same
  projections as before.
- **Defaults stay fixed.** `geometry_default` / `defaultGeometry` and the
  default values of optional fields (such as `accuracy`, `COR`,
  `rotDetector` or `mode`) do not change.
- **Validation may get stricter.** A new check may reject a geometry that
  was always inconsistent (for example `sVoxel` not equal to
  `nVoxel * dVoxel`), but only with an error message that says what to fix.

A change to a convention (for example the direction of rotation or the
position of the source at angle 0) is only made in a major release. Before
that, the old behaviour is deprecated and the documentation explains how to
convert an existing geometry.

## Numerical results

TIGRE is used to compare methods, so changes to results are made openly:

- A **bug fix** may change results, even in a patch release. The changelog
  says what changed and why (for example, the Parker and Wang weighting
  changes in 3.2.0).
- A **default value** of an algorithm option (such as a regularisation
  parameter or a step size) may change in a minor release, when the new value
  is clearly better. It is listed under *Changed* in the changelog. To keep
  exactly the previous behaviour, pass the option explicitly.
- Results are not bit-for-bit reproducible across GPUs, CUDA versions or
  numbers of GPUs: floating-point operations on the GPU are not always done
  in the same order. Differences are at the level of rounding errors.

## How deprecation works

When a function, an option or a behaviour is going to be removed or changed:

1. It is marked as deprecated in a minor release: the documentation says so
   and what to use instead, and the changelog lists it under *Deprecated*.
2. Using it gives a warning: a `FutureWarning` in Python (always shown), or
   a MATLAB warning with identifier
   `TIGRE:deprecated`, which you can silence with
   `warning('off', 'TIGRE:deprecated')`.
3. It keeps working for at least one more minor release and at least six
   months, and is removed in the next major release after that.

## Supported environments

- **Python:** the versions that are still supported upstream. A Python
  version is dropped in a minor release once it reaches its end of life (as
  Python 3.10 was in 3.2.0).
- **CUDA, compilers and MATLAB:** see [](../install/index.md) for the versions
  currently tested. Support for old CUDA versions or GPU architectures may be
  dropped in a minor release when it stops TIGRE from building with current
  tools; the changelog says when.
