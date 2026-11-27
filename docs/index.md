---
html_theme.sidebar_secondary.remove: true
og:description: "TIGRE is a GPU toolbox for 3D tomographic (CT) reconstruction in MATLAB and Python: any scanning geometry, multi-GPU projectors and a wide range of iterative algorithms."
myst:
  html_meta:
    google-site-verification: "Qcj7Q89d0NZ0WwCxvdHUY1TSH7blFCSGS78LsGrpN20"
    description: "TIGRE is a GPU toolbox for 3D tomographic (CT) reconstruction in MATLAB and Python: any scanning geometry, multi-GPU projectors and a wide range of iterative algorithms."
---

# TIGRE

<div class="tigre-band">
<div class="tigre-band-inner">
<div class="tigre-band-head">

```{image} _static/logo/tigre-logo-dark.svg
:alt: TIGRE
:class: tigre-hero-logo
```

<p class="tigre-tagline">Fast, flexible 3D CT reconstruction on GPUs, for MATLAB and Python.</p>

</div>
<div class="tigre-band-body">

TIGRE, the Tomographic Iterative GPU-based Reconstruction toolbox,
reconstructs images from any scanning geometry with multi-GPU projectors
and a wide range of iterative algorithms.

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

<p class="tigre-news">

**New in {{ latest_version }}** ({{ latest_date }}): faster projectors and
backprojectors on one and several GPUs. Read the
[release notes](changelog.md) or {{ latest_download }}.

</p>

</div>

<figure class="tigre-compare" aria-label="Reconstructions compared">
<div class="tigre-compare-tabs" role="tablist" aria-label="Choose a scan">
<button type="button" role="tab" id="tab-micro" aria-controls="panel-micro" aria-selected="true">Micro-CT</button>
<button type="button" role="tab" id="tab-synchrotron" aria-controls="panel-synchrotron" aria-selected="false" tabindex="-1">Synchrotron</button>
<button type="button" role="tab" id="tab-clinical" aria-controls="panel-clinical" aria-selected="false" tabindex="-1">Clinical CBCT</button>
</div>
<div class="tigre-compare-panel" role="tabpanel" id="panel-micro" aria-labelledby="tab-micro">
<div class="tigre-compare-frame" style="--pos: 82%">
<img class="dark-light" src="_static/hero/coffee-cgls.webp" width="1000" height="1000" alt="The same slice reconstructed with CGLS: the same structure, with less noise">
<img class="tigre-compare-left-img dark-light" src="_static/hero/coffee-fdk.webp" width="1000" height="1000" alt="A slice through a coffee bean reconstructed with FDK: fine structure with grainy noise">
<div class="tigre-compare-zoom" aria-hidden="true">
<img src="_static/hero/coffee-zoom-cgls.webp" width="744" height="744" alt="" loading="lazy">
<img class="tigre-compare-zoom-left" src="_static/hero/coffee-zoom-fdk.webp" width="744" height="744" alt="" loading="lazy">
</div>
<span class="tigre-compare-name tigre-compare-left" aria-hidden="true">FDK</span>
<span class="tigre-compare-name tigre-compare-right" aria-hidden="true">CGLS</span>
<input class="tigre-compare-range" type="range" min="0" max="100" value="82" step="1" aria-label="Move the divider between FDK (left) and CGLS (right)">
<span class="tigre-compare-handle" aria-hidden="true"></span>
</div>
<p class="tigre-compare-caption">A coffee bean: 3340 × 3340 × 900 voxels from 2134 projections, larger than GPU memory. FDK against 15 iterations of CGLS. The box at the top right magnifies the marked area, and grows while you drag.</p>
</div>
<div class="tigre-compare-panel" role="tabpanel" id="panel-synchrotron" aria-labelledby="tab-synchrotron" hidden>
<div class="tigre-compare-frame" style="--pos: 50%">
<img class="dark-light" src="_static/hero/synchrotron-ospocs.webp" width="490" height="490" alt="The same slice reconstructed with OS-ASD-POCS: clear, sharp grains" loading="lazy">
<img class="tigre-compare-left-img dark-light" src="_static/hero/synchrotron-fbp.webp" width="490" height="490" alt="A synchrotron slice reconstructed with FBP: almost pure noise" loading="lazy">
<span class="tigre-compare-name tigre-compare-left" aria-hidden="true">FBP</span>
<span class="tigre-compare-name tigre-compare-right" aria-hidden="true">OS-ASD-POCS</span>
<input class="tigre-compare-range" type="range" min="0" max="100" value="50" step="1" aria-label="Move the divider between FBP (left) and OS-ASD-POCS (right)">
<span class="tigre-compare-handle" aria-hidden="true"></span>
</div>
<p class="tigre-compare-caption">A synchrotron scan with the fast, low-dose settings of an in-situ experiment. FBP against 20 iterations of OS-ASD-POCS. From <a href="https://doi.org/10.1088/2631-8695/adbb3a">Biguri et al. (2025)</a>, CC BY 4.0.</p>
</div>
<div class="tigre-compare-panel" role="tabpanel" id="panel-clinical" aria-labelledby="tab-clinical" hidden>
<div class="tigre-compare-frame" style="--pos: 50%; aspect-ratio: 1436 / 752">
<img class="dark-light" src="_static/hero/lung-awpcsd.webp" width="1436" height="752" alt="A coronal chest slice from a low-dose cone-beam scan reconstructed with AwPCSD: lungs, vessels and a tumour are clear, with little noise" loading="lazy">
<img class="tigre-compare-left-img dark-light" src="_static/hero/lung-fdk.webp" width="1436" height="752" alt="The same slice reconstructed with FDK: the same anatomy with strong grainy noise" loading="lazy">
<span class="tigre-compare-name tigre-compare-left" aria-hidden="true">FDK</span>
<span class="tigre-compare-name tigre-compare-right" aria-hidden="true">AwPCSD</span>
<input class="tigre-compare-range" type="range" min="0" max="100" value="50" step="1" aria-label="Move the divider between FDK (left) and AwPCSD (right)">
<span class="tigre-compare-handle" aria-hidden="true"></span>
</div>
<p class="tigre-compare-caption">Low-dose cone-beam CT of a lung cancer patient for image-guided radiotherapy. FDK against AwPCSD, an edge-preserving iterative method. From <a href="https://doi.org/10.1088/1361-6560/aa7675">Biguri et al. (2017)</a>, Physics in Medicine &amp; Biology.</p>
</div>
<figcaption class="tigre-compare-hint">Each method is named on its side of the orange divider. Drag it, or use the arrow keys, to compare.</figcaption>
</figure>

</div>
</div>

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
:gutter: 4
:class-container: tigre-why

:::{grid-item}
**[Any geometry](user-guide/geometry.md)**

Cone beam, parallel beam, helical, tomosynthesis, C-arm and curved
detectors. Geometric parameters are defined per projection, not per scan.
:::

:::{grid-item}
**[State-of-the-art speed](technical-notes/speed-3.3.1/index.md)**

Optimised CUDA projectors and backprojectors, among the fastest available,
that split problems larger than GPU memory across one or several GPUs.
:::

:::{grid-item}
**[Many algorithms](user-guide/algorithms.md)**

FDK/FBP, the SART family, Krylov methods (CGLS, LSQR, LSMR, GMRES, hybrid),
MLEM, total-variation and prior-image methods, and more.
:::

:::{grid-item}
**[MATLAB and Python](user-guide/interoperability.md)**

The same toolbox and the same demos in both languages, and PyTorch
operators with gradients for learned reconstruction.
:::

:::{grid-item}
**[Real data](user-guide/loading-data.md)**

Loaders for Nikon, Bruker, YXLON, Diondo and Varian scanners, DXchange
synchrotron data and DICOM (MATLAB).
:::

:::{grid-item}
**[Widely used](papers.md)**

Used in {{ papers_count }} published papers, in medicine, synchrotron and
industrial imaging, particle therapy and reconstruction research.
:::
::::

## From real scans

<div class="tigre-strip">
<a href="gallery/index.html#a-clinical-scanner-varian-edge"><img src="_static/home/clinical-cbct.webp" width="640" height="480" alt="Head scanned on a clinical Varian cone-beam CT" loading="lazy"><span>Clinical cone-beam CT</span></a>
<a href="gallery/index.html#synchrotron-in-situ-settings"><img src="_static/home/synchrotron.webp" width="640" height="480" alt="Synchrotron scan reconstructed with several algorithms" loading="lazy"><span>Synchrotron</span></a>
<a href="gallery/index.html#micro-ct-a-helical-scan"><img src="_static/home/helical-micro-ct.webp" width="640" height="480" alt="Helical micro-CT scan of a battery" loading="lazy"><span>Helical micro-CT</span></a>
<a href="gallery/index.html#neutron-tomography"><img src="_static/home/neutron.webp" width="640" height="480" alt="Neutron tomography of a padlock" loading="lazy"><span>Neutron tomography</span></a>
</div>

[More in the Gallery](gallery/index.md)

## Get started

::::{grid} 1 1 3 3
:gutter: 3
:class-container: tigre-start

:::{grid-item-card} Install
:link: install/index
:link-type: doc
Prebuilt wheels for Python on Linux and Windows, or build TIGRE for MATLAB
and Python from source.
:::

:::{grid-item-card} Learn
:link: user-guide/index
:link-type: doc
The user guide: geometry, algorithms, real data, and what to check when a
reconstruction looks wrong.
:::

:::{grid-item-card} Reference
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
