# Gallery

Reconstructions made with TIGRE, mostly from real scans, across medical,
synchrotron, micro-CT, neutron and proton imaging. For papers that use TIGRE
in each field, see [](../papers.md).

## Medical cone-beam CT

A RANDO head phantom, which represents an adult head, scanned at The
Christie Hospital in Manchester with 360 equally spaced projections over a
full rotation. The data have mechanical offsets and high noise, from head
CBCT settings with low X-ray intensity.

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/randofull.png
:alt: RANDO head phantom reconstructed with four algorithms from 360 projections

RANDO head from 360 projections. Left to right, top to bottom: FDK, CGLS,
OS-SART and OS-ASD-POCS. Display window [0, 0.05].
```

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/rando90.png
:alt: RANDO head phantom reconstructed with four algorithms from 90 projections

The same scan with only 90 projections, same algorithms and display window.
```

Thanks to The Christie Hospital for sharing this dataset. Source: Ander
Biguri's [PhD thesis](https://github.com/AnderBiguri/PhDThesis), page 93.

### A clinical scanner: Varian Edge

A head scanned over a 200° arc on a Varian Edge machine at Beijing Cancer
Hospital.

```{figure} images/v3-varian-head.jpg
:alt: Head slices reconstructed with FDK, OS-SART and OS-ASD-POCS

Rows are different slices; columns are FDK, OS-SART (50 iterations) and
OS-ASD-POCS (50 iterations). Display window [0, 0.05].
```

Thanks to Yi Liu for the data. Source: [Biguri et al. (2025), TIGRE v3](https://doi.org/10.1088/2631-8695/adbb3a), figure 2 ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

### Low-dose lung imaging

Image-guided radiotherapy for lung cancer images the patient repeatedly, so
every scan should use as little dose as possible. Regularised iterative
algorithms remove much of the noise that FDK leaves in low-dose cone-beam CT.

```{figure} images/lung-cbct.jpg
:alt: Coronal lung slices: reference image, FDK, OS-SART and AwPCSD

Reference image, FDK, OS-SART and AwPCSD.
```

Source: [Biguri et al. (2017)](https://doi.org/10.1088/1361-6560/aa7675),
Physics in Medicine & Biology.

## Cryo soft X-ray tomography

Trypanosomes (single-celled parasites that cause sleeping sickness), frozen
and imaged with soft X-rays at beamline B24 of the Diamond Light Source.
Soft X-rays reach nanometre resolution, but the scans are very noisy and
cover less than 120 degrees.

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TV.png
:alt: Trypanosoma reconstructed with FBP, OS-SART and OS-ASD-POCS
```

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TVz1.png
:alt: Zoomed region of the Trypanosoma reconstructions
```

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TVz2.png
:alt: Second zoomed region of the Trypanosoma reconstructions

Columns: FBP, OS-SART (20 iterations) and OS-ASD-POCS (20 iterations, 20 TV
iterations each). The red squares in the first image mark the zoomed areas.
```

Thanks to the Diamond Light Source, especially Michele Darrow and Imanol
Luengo, for sharing this dataset. Source: PhD thesis, page 101.

## Synchrotron: in-situ settings

A synchrotron scan acquired with the fast, low-dose settings of a dynamic
in-situ experiment, where every algorithm has to work with noisy data.

```{figure} images/v3-synchrotron.jpg
:alt: Synchrotron slice reconstructed with six algorithms

Left to right, top to bottom: FBP (17 s), LSQR with 20
iterations (249 s), FISTA with 80 iterations (1330 s), FBP with a Hamming
filter (17 s), SIRT with 150 iterations (1330 s) and OS-ASD-POCS with 20
iterations (856 s).
```

Source: [Biguri et al. (2025), TIGRE v3](https://doi.org/10.1088/2631-8695/adbb3a), figure 3 ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

## Micro-CT: SophiaBeads

The [SophiaBeads dataset](https://doi.org/10.5281/zenodo.16474) was acquired
for testing reconstruction methods: a 25 mm plastic tube filled with soda-lime
glass beads 2.5 mm across. Here it is reconstructed from 256 projections with
FDK and CGLS (15 iterations).

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/sophiaFDKCGLS.png
:alt: SophiaBeads reconstructed with FDK and CGLS

FDK and CGLS, 256 projections. The red line marks the profile below. Display
window [0, 0.15].
```

```{figure} https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/spophiaprofile.png
:alt: Intensity profile along the red line for FDK and CGLS
```

Source: PhD thesis, page 96.

## Micro-CT: a very large reconstruction

A coffee bean scanned with a shifted detector panel, reconstructed at
3340 × 3340 × 900 voxels from 2134 projections of 3780 × 900 pixels. It shows
that TIGRE can reconstruct images much larger than GPU memory, even with
CGLS, which needs a lot of memory.

```{figure} ../../Frontispiece/coffee.png
:alt: Coffee bean reconstructed with FDK and CGLS

Coffee bean with [FDK](https://github.com/CERN/TIGRE/blob/master/Frontispiece/fdk.png?raw=true)
and [CGLS](https://github.com/CERN/TIGRE/blob/master/Frontispiece/cgls.png?raw=true)
(15 iterations). The image is heavily compressed; follow the links for the
full-resolution originals.
```

Source: [Biguri et al. (2020)](https://doi.org/10.1016/j.jpdc.2020.07.004),
figure 10.

## Micro-CT: a helical scan

An Energizer AAAA battery scanned in helical mode on a diondo d5 scanner at
the µ-VIS X-ray Imaging Centre, University of Southampton.

```{figure} images/v3-battery.jpg
:alt: Battery reconstructed with FDK and with FISTA-TV

Top: FDK. Bottom: FISTA with TV regularisation (20 iterations).
```

Source: [Biguri et al. (2025), TIGRE v3](https://doi.org/10.1088/2631-8695/adbb3a), figure 6 ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

## Neutron tomography

A padlock imaged with neutrons at the Thailand Institute of Nuclear
Technology (TINT).

```{figure} images/v3-neutron-padlock.jpg
:alt: Padlock reconstructed with FDK, SART, SART-TV and SIRT

Left to right: FDK with a Shepp-Logan filter, SART, SART-TV and SIRT, the
iterative ones with 100 iterations each.
```

Source: [Biguri et al. (2025), TIGRE v3](https://doi.org/10.1088/2631-8695/adbb3a), figure 5 ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

## Proton CT

Simulated proton CT of the high-resolution module of the Catphan phantom,
reconstructed with TIGRE's MATLAB proton CT extension (see
[](../user-guide/proton-ct.md)).

```{figure} images/v3-proton-ct.jpg
:alt: Catphan high-resolution module reconstructed with FBP, ASD-POCS and OSEM

Left to right: FBP, ASD-POCS and OSEM.
```

Source: [Biguri et al. (2025), TIGRE v3](https://doi.org/10.1088/2631-8695/adbb3a), figure 4 ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

## Designing scan trajectories

TIGRE's geometry can be set separately for every projection, so it can
simulate scans along any source and detector path. Small dental CBCT units
have a small field of view; this work uses TIGRE to design arc-based
trajectories that move the scan towards chosen regions outside the standard
field of view, such as the temporal bones or the spine, while using fewer
projections.

```{figure} images/trajectory.png
:alt: Two arc-based dental CBCT trajectories (top) and the corresponding head reconstructions (bottom)

Top: trajectories designed for two regions of interest (masks 2 and 3).
Bottom: reconstructions of a head from scans along them.
```

Source: [Islam et al. (2025)](https://doi.org/10.1016/j.csbj.2025.06.013),
Computational and Structural Biotechnology Journal
([CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/)).
