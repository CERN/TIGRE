# Reconstructing a real scan, step by step

This tutorial reconstructs a real micro-CT scan from start to finish: download
the data, load it, check it, fix the centre of rotation, and reconstruct the
whole object. It uses Python; the MATLAB functions have the same names
([](loading-data.md)).

The data are the **SophiaBeads** datasets ([Coban and McDonald,
2015](https://doi.org/10.5281/zenodo.16474)): a plastic tube, 25 mm across,
filled with glass beads about 2.5 mm across, scanned on a Nikon micro-CT
scanner at the University of Manchester. The same sample was scanned with
2048, 1024, ... down to 64 projections. This tutorial uses the 64-projection
scan, the smallest.

```{admonition} Download size
:class: warning

The 64-projection dataset is a **0.47 GB** zip file, 0.5 GB unzipped. You
need about 1 GB of free disk space while it unzips, about 6 GB of RAM, and an
NVIDIA GPU. The larger datasets on the same Zenodo record are much bigger: up
to 16 GB for 2048 projections.

The data are not part of TIGRE. They are licensed CC BY-SA 4.0; if you publish
work that uses them, cite the dataset
([doi:10.5281/zenodo.16474](https://doi.org/10.5281/zenodo.16474)) as the
Zenodo record asks.
```

## 1. Download the data

Download `SophiaBeads_64_averaged.zip` from the
[Zenodo record](https://zenodo.org/records/16474) and unzip it. From the
command line, the [`zenodo-get`](https://pypi.org/project/zenodo-get/) package
does the same:

```bash
pip install zenodo-get
zenodo_get 16474 -g "*_64_*"
```

The unzipped folder, `SophiaBeads_64_averaged`, holds the projections (TIFF
images), the angles (`.ang`) and the scanner settings (`.xtekct`).

The script `docs/tools/make_sophiabeads_tutorial.py` in the TIGRE repository
downloads the data and runs every step of this page; it made the pictures
below.

## 2. Load it

```python
import numpy as np
import tigre
import tigre.algorithms as algs
from tigre.utilities.io import NikonDataLoader

proj, geo, angles = NikonDataLoader("SophiaBeads_64_averaged")
print(proj.shape, proj.min(), proj.max())
print(geo)
```

```{admonition} TIGRE 3.3.1 and earlier
:class: note

These versions fail on this dataset with `TypeError: must be real number, not
str`, because they read the `.ang` file's angles as text. It is fixed in the
next version; until then, the function `load_nikon` in the tutorial script
works around it.
```

The loader reads the scanner settings and returns three things:

- `proj`: 63 projections of 2000 × 2000 pixels, already converted to line
  integrals ($-\log(I/I_0)$). That is 1 GB in memory. (The scanner settings
  say 63 projections; the folder has one more image, which is not used.)
- `geo`: the geometry. The source is 80.6 mm from the rotation axis
  (`geo.DSO`) and 1007 mm from the detector (`geo.DSD`), so the object is
  magnified 12.5 times. The detector pixels are 0.2 mm (`geo.dDetector`), so
  the voxels are 0.2 / 12.5 = 0.016 mm (`geo.dVoxel`). The scanner's
  suggested image is 2000 × 2000 × 2000 voxels (`geo.nVoxel`).
- `angles`: 63 angles over a full turn, in radians.

Look at a projection before reconstructing anything:

```python
tigre.plotproj(proj, angles)
```

```{image} images/sophiabeads/projection.png
:alt: One projection of the bead tube: a bright vertical band of overlapping circles
:width: 360px
```

Air is close to 0 and the beads are bright, as line integrals should be. If
the picture were the other way round, the data would still be intensities
([](artefacts.md#log-and-flat-field)).

## 3. Check the geometry on a few slices

The full image is 2000³ voxels: 32 GB in single precision. Before
reconstructing all of it, reconstruct a few slices in the middle. That takes
about a second and shows whether the geometry is right.

```python
import copy

test = copy.deepcopy(geo)
test.nVoxel = np.array([8, geo.nVoxel[1], geo.nVoxel[2]])
test.sVoxel = test.nVoxel * test.dVoxel
img = algs.fdk(proj, test, angles)
tigre.plotimg(img, dim="z")
```

The left picture below is the result. The beads have soft, doubled edges and
dark arcs between them. That is the signature of a wrong **centre of
rotation**: the rotation axis does not project exactly onto the centre of the
detector ([](artefacts.md#centre-of-rotation)).

## 4. Fix the centre of rotation

`computeCOR` estimates the shift from the data, by finding where the rays
measured half a turn apart agree:

```python
geo.COR = tigre.computeCOR(proj, geo, angles)
test.COR = geo.COR
img = algs.fdk(proj, test, angles)
```

Here it finds `geo.COR = 0.234` mm, about 15 voxels. With it, the beads are
round and sharp (right):

```{image} images/sophiabeads/cor.png
:alt: The same slice before and after setting the centre of rotation; after it, the beads are sharp
```

```{image} images/sophiabeads/cor_zoom.png
:alt: Close-up of the same slices: doubled bead edges before, clean round beads after
```

Many micro-CT scans need this correction. Check it whenever the scanner's file
does not give the centre of rotation (this one gives 0).

## 5. Reconstruct the whole object

The full 2000³ image does not fit on most GPUs. TIGRE can split it between
several GPUs and reconstruct it in pieces ([](gpus.md)), but for a first look
a quarter of the resolution is enough: 500³ voxels of 0.064 mm, the same field
of view, 0.5 GB.

```python
small = copy.deepcopy(geo)
small.nVoxel = (geo.nVoxel // 4).astype(int)
small.dVoxel = geo.sVoxel / small.nVoxel

img_fdk = algs.fdk(proj, small, angles)
img_cgls = algs.cgls(proj, small, angles, 15)
```

```{image} images/sophiabeads/full_axial.png
:alt: Axial slice of the whole tube reconstructed with FDK and with CGLS
```

```{image} images/sophiabeads/full_coronal.png
:alt: Vertical slice through the tube reconstructed with FDK and with CGLS
```

With only 63 projections, FDK (left) leaves streaks between and around the
beads ([](artefacts.md#too-few-projections)). CGLS (right) has much weaker
streaks and a smoother background, at the cost of slightly softer bead edges.
More iterations sharpen it, but with noisy data they eventually add noise too
([](artefacts.md#number-of-iterations)); [](algorithms.md) compares the
algorithms. The bright arcs at the top and bottom of the vertical CGLS slice
are the edges of the cone of rays: those corners of the image are seen by only
some of the projections.

On one desktop GPU, FDK takes about a second here and 15 iterations of CGLS
about a minute.

## 6. Save the result

The images are NumPy arrays, indexed `[z, y, x]`. Save them in any format,
for example:

```python
np.save("sophiabeads_cgls.npy", img_cgls)
```

or as a stack of TIFF images with `tifffile.imwrite("sophiabeads_cgls.tif",
img_cgls)`, which ImageJ/Fiji and most other tools can open.

## Next steps

- Try a larger SophiaBeads dataset (more projections): the streaks in FDK
  fade as the number of projections grows.
- [](geometry-walkthrough.md) explains each geometry check in more detail, for
  scanners without a loader.
- [](artefacts.md) shows the other common mistakes and their fixes.
