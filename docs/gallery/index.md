# Gallery

Reconstructions of real scans made with TIGRE. They show TIGRE's algorithms on real scans from different
kinds of CT. For papers that use TIGRE in each field, see
[](../papers.md).

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
