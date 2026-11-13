---
html_theme.sidebar_secondary.remove: true
---

# TIGRE

::::{grid} 1 1 2 2
:gutter: 4
:class-container: tigre-hero

:::{grid-item}
:columns: 12 12 6 6

<p class="tigre-tagline">The Tomographic Iterative GPU-based Reconstruction toolbox: fast, flexible 3D CT reconstruction for MATLAB and Python.</p>

TIGRE reconstructs tomographic images from any scanning geometry, using GPU
implementations (including multi-GPU) of the forward and backprojection
operators and a wide range of iterative algorithms.

```{button-ref} install/index
:ref-type: doc
:color: primary
:class: sd-mr-2

Install TIGRE
```

```{button-link} https://github.com/CERN/TIGRE
:color: secondary
:outline:

View on GitHub
```
:::

:::{grid-item}
:columns: 12 12 6 6

```{figure} ../Frontispiece/coffee.png
:alt: Two reconstructions of a coffee bean slice, side by side

Micro-CT of a coffee bean, 3340 × 3340 × 900 voxels from 2134 projections,
reconstructed with FDK and with CGLS.
```
:::
::::

## Why TIGRE

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} Any geometry
Cone beam, parallel beam, helical, tomosynthesis, C-arm and curved
detectors. Geometric parameters are defined per projection, not per scan.
:::

:::{grid-item-card} Fast on GPUs
CUDA projectors and backprojectors that split large problems across
multiple GPUs automatically.
:::

:::{grid-item-card} Many algorithms
FDK/FBP, SART-family, Krylov subspace methods (CGLS, LSQR, LSMR, GMRES),
MLEM, total-variation regularised methods and more.
:::

:::{grid-item-card} MATLAB and Python
The same toolbox and the same demos in both languages.
:::

:::{grid-item-card} Real data
Loaders for Nikon, Varian, Philips (DICOM) and other scanner formats.
:::

:::{grid-item-card} Open source
BSD-licensed, developed at the University of Bath and CERN with
contributions from the community.
:::
::::

## Get started

::::{grid} 1 2 2 4
:gutter: 3

:::{grid-item-card} Install
:link: install/index
:link-type: doc
MATLAB and Python installation.
:::

:::{grid-item-card} User guide
:link: user-guide/index
:link-type: doc
Geometry, projectors, algorithms and more.
:::

:::{grid-item-card} Gallery
:link: gallery/index
:link-type: doc
Worked examples from the demos.
:::

:::{grid-item-card} Cite TIGRE
:link: citing
:link-type: doc
How to cite TIGRE in your work.
:::
::::

```{toctree}
:hidden:
:maxdepth: 2

install/index
user-guide/index
gallery/index
api/index
citing
papers
```
