---
html_theme.sidebar_secondary.remove: true
og:description: "TIGRE is a GPU toolbox for 3D tomographic (CT) reconstruction in MATLAB and Python: any scanning geometry, multi-GPU projectors and a wide range of iterative algorithms."
myst:
  html_meta:
    description: "TIGRE is a GPU toolbox for 3D tomographic (CT) reconstruction in MATLAB and Python: any scanning geometry, multi-GPU projectors and a wide range of iterative algorithms."
---

# TIGRE

::::{grid} 1 1 2 2
:gutter: 4
:class-container: tigre-hero

:::{grid-item}
:columns: 12 12 6 6

```{image} _static/logo/tigre-logo.svg
:alt: TIGRE
:class: tigre-hero-logo only-light
```

```{image} _static/logo/tigre-logo-dark.svg
:alt: TIGRE
:class: tigre-hero-logo only-dark
```

<p class="tigre-tagline">The Tomographic Iterative GPU-based Reconstruction toolbox: fast, flexible 3D CT reconstruction for MATLAB and Python.</p>

TIGRE reconstructs tomographic images from any scanning geometry, using GPU
implementations (including multi-GPU) of the forward and backprojection
operators and a wide range of iterative algorithms.

<div class="tigre-actions">

```{button-ref} install/index
:ref-type: doc
:color: primary

Install TIGRE
```

```{button-ref} demos/index
:ref-type: doc
:color: secondary
:outline:

See the demos
```

```{button-link} https://github.com/CERN/TIGRE
:color: secondary
:outline:

View on GitHub
```

</div>

<div class="tigre-news">

**New: TIGRE {{ latest_version }}** ({{ latest_date }}). Faster projectors and
backprojectors on one and several GPUs.
[Release notes](changelog.md) ·
{{ latest_download }}

</div>
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

## A reconstruction in a few lines

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```python
import numpy as np
import tigre
import tigre.algorithms as algs
from tigre.utilities import sample_loader

geo = tigre.geometry_default(high_resolution=False)   # a cone-beam scanner
angles = np.linspace(0, 2 * np.pi, 100)
head = sample_loader.load_head_phantom(geo.nVoxel)

proj = tigre.Ax(head, geo, angles)                     # simulate the scan
img_fdk = algs.fdk(proj, geo, angles)                  # filtered backprojection
img_ossart = algs.ossart(proj, geo, angles, 20)        # 20 iterations of OS-SART
```
:::

:::{tab-item} MATLAB
:sync: matlab

```matlab
geo = defaultGeometry('nVoxel', [128; 128; 128]);     % a cone-beam scanner
angles = linspace(0, 2*pi, 100);
head = headPhantom(geo.nVoxel);

proj = Ax(head, geo, angles);                         % simulate the scan
imgFDK = FDK(proj, geo, angles);                      % filtered backprojection
imgOSSART = OS_SART(proj, geo, angles, 20);           % 20 iterations of OS-SART
```
:::
::::

Swap the phantom for your own projections, or use one of the
[scanner loaders](user-guide/loading-data.md), and change the geometry to
match your system. The [demos](demos/index.rst) go through every step.

## Why TIGRE

::::{grid} 1 2 2 3
:gutter: 3
:class-container: tigre-features

:::{grid-item-card} {octicon}`sync;1.4em;tigre-icon` Any geometry
:link: user-guide/geometry
:link-type: doc
Cone beam, parallel beam, helical, tomosynthesis, C-arm and curved
detectors. Geometric parameters are defined per projection, not per scan.
:::

:::{grid-item-card} {octicon}`zap;1.4em;tigre-icon` State-of-the-art speed
:link: technical-notes/speed-3.3.1/index
:link-type: doc
Optimised CUDA projectors and backprojectors, among the fastest available,
that split problems larger than GPU memory across one or several GPUs.
:::

:::{grid-item-card} {octicon}`stack;1.4em;tigre-icon` Many algorithms
:link: user-guide/algorithms
:link-type: doc
FDK/FBP, the SART family, Krylov methods (CGLS, LSQR, LSMR, GMRES, hybrid),
MLEM, total-variation and prior-image methods, and more.
:::

:::{grid-item-card} {octicon}`code;1.4em;tigre-icon` MATLAB and Python
:link: user-guide/interoperability
:link-type: doc
The same toolbox and the same demos in both languages, and PyTorch
operators with gradients for learned reconstruction.
:::

:::{grid-item-card} {octicon}`database;1.4em;tigre-icon` Real data
:link: user-guide/loading-data
:link-type: doc
Loaders for Nikon, Bruker, YXLON, Diondo and Varian scanners, DXchange
synchrotron data and DICOM (MATLAB).
:::

:::{grid-item-card} {octicon}`people;1.4em;tigre-icon` Widely used
:link: papers
:link-type: doc
Used in {{ papers_count }} published papers, in medicine, synchrotron and
industrial imaging, particle therapy and reconstruction research.
:::
::::

## From real scans

::::{grid} 2 2 4 4
:gutter: 2
:class-container: tigre-strip

:::{grid-item-card}
:link: gallery/index
:link-type: doc
```{image} gallery/images/v3-varian-head.jpg
:alt: Head scanned on a clinical Varian cone-beam CT
```

Clinical cone-beam CT
:::

:::{grid-item-card}
:link: gallery/index
:link-type: doc
```{image} gallery/images/v3-synchrotron.jpg
:alt: Synchrotron scan reconstructed with several algorithms
```

Synchrotron
:::

:::{grid-item-card}
:link: gallery/index
:link-type: doc
```{image} gallery/images/v3-battery.jpg
:alt: Helical micro-CT scan of a battery
```

Helical micro-CT
:::

:::{grid-item-card}
:link: gallery/index
:link-type: doc
```{image} gallery/images/v3-neutron-padlock.jpg
:alt: Neutron tomography of a padlock
```

Neutron tomography
:::
::::

[More in the Gallery](gallery/index.md)

## Get started

::::{grid} 1 1 3 3
:gutter: 3
:class-container: tigre-features tigre-start

:::{grid-item-card} {octicon}`download;1.4em;tigre-icon` Install
:link: install/index
:link-type: doc
Prebuilt wheels for Python on Linux, or build TIGRE for MATLAB and Python
from source.
:::

:::{grid-item-card} {octicon}`mortar-board;1.4em;tigre-icon` Learn
:link: user-guide/index
:link-type: doc
The user guide: geometry, algorithms, real data, and what to check when a
reconstruction looks wrong.
:::

:::{grid-item-card} {octicon}`book;1.4em;tigre-icon` Reference
:link: api/index
:link-type: doc
Every function and option, in Python and MATLAB.
:::
::::

## Who makes TIGRE

TIGRE is open source under the
[BSD 3-clause licence](https://github.com/CERN/TIGRE/blob/master/LICENSE.txt).
The copyright and the licence are held by the University of Bath and CERN,
where TIGRE was first developed. Ander Biguri has led its development since,
at the University of Southampton and now at the University of Cambridge, with
[contributions from many people](https://github.com/CERN/TIGRE#contributors-).
The Universities of Southampton and Cambridge do not hold the copyright or
the licence.

TIGRE has been funded by CERN, the EPSRC, the European Metrology Research
Programme and the Accelerate Programme for Scientific Discovery; see the
[acknowledgements](acknowledgements.md).

```{toctree}
:hidden:
:maxdepth: 2

install/index
user-guide/index
demos/index
gallery/index
api/index
citing
papers
technical-notes/index
ecosystem
roadmap
help
contributing
changelog
acknowledgements
```
