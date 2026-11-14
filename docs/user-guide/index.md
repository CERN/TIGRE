# User guide

This guide explains how TIGRE works and how to use it on your own data. It
assumes TIGRE is installed ([](../install/index.md)).

## A first reconstruction

A TIGRE reconstruction needs three things: a **geometry** describing the
scanner, the **angles** of the projections, and the **projections**
themselves. This example simulates projections of a head phantom and
reconstructs them two ways.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python
```python
import numpy as np
import tigre
import tigre.algorithms as algs
from tigre.utilities import CTnoise, sample_loader

# 1. Geometry and angles
geo = tigre.geometry_default(high_resolution=False)
angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)

# 2. Projections: here simulated from a phantom, with realistic noise
head = sample_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(head, geo, angles)
proj = CTnoise.add(proj, Poisson=1e5, Gaussian=np.array([0, 10]))

# 3. Reconstruct
img_fdk = algs.fdk(proj, geo, angles)
img_sart = algs.ossart(proj, geo, angles, 20)

tigre.plotimg(np.concatenate([img_fdk, img_sart], axis=1), dim="z")
```
:::

:::{tab-item} MATLAB
:sync: matlab
```matlab
% 1. Geometry and angles
geo = defaultGeometry('nVoxel', [128; 128; 128]);
angles = linspace(0, 2*pi, 100);

% 2. Projections: here simulated from a phantom, with realistic noise
head = headPhantom(geo.nVoxel);
proj = Ax(head, geo, angles);
proj = addCTnoise(proj);

% 3. Reconstruct
img_fdk = FDK(proj, geo, angles);
img_sart = OS_SART(proj, geo, angles, 20);

plotImg([img_fdk img_sart], 'Dim', 'Z');
```
:::
::::

With real data, step 2 is loading your scan instead ([](loading-data.md)).

## Guide

[](geometry.md)
: The geometry fields, the coordinate system, and how to describe offsets,
  helical scans and other scanners. Start here.

[](projectors.md)
: The forward and backprojection operators everything else is built on.

[](algorithms.md)
: Which algorithm to use, and the options they share.

[](loading-data.md)
: Reading data from Nikon, Bruker, Varian and other scanners.

[](gpus.md)
: Choosing GPUs, how large problems are split, and memory use.

[](troubleshooting.md)
: Common errors and how to fix them.

[](performance-tuning.md)
: Source-level changes for a little more speed.

[](proton-ct.md)
: The MATLAB extension for proton radiography.

[](further-reading.md)
: Papers and the thesis behind TIGRE.

Every feature also has a demo; the [demos](../demos/index.rst) show
them all, in Python and MATLAB.

```{toctree}
:hidden:

geometry
projectors
algorithms
loading-data
gpus
troubleshooting
performance-tuning
proton-ct
further-reading
```
