---
og:description: "Reference for TIGRE's MATLAB functions: Ax and Atb, geometry, reconstruction algorithms, scanner data loaders and utilities."
myst:
  html_meta:
    description: "Reference for TIGRE's MATLAB functions: Ax and Atb, geometry, reconstruction algorithms, scanner data loaders and utilities."
---

# MATLAB API

Generated from each function's help text, which is what `help NAME` shows in
MATLAB. Arrays are `(x, y, z)` for images and `(u, v, angle)` for projections,
and options are given as name-value pairs, for example
`OS_SART(proj, geo, angles, 20, 'BlockSize', 10)`.

```{toctree}
:maxdepth: 1

Projection operators <operators>
Geometry <geometry>
Algorithms <algorithms>
Utilities <utilities>
Scanner data loaders <io>
```
