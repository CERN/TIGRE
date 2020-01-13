Gallery: TIGRE algorithms in different CT applications. 
=======

## Conventional medical CBCT

**Source of image:** Ander Biguri's PhD thesis, page 93.

**Description:** This dataset is a RANDO head phantom that represents a head of an adult. The
dataset was acquired in The Christie Hospital in Manchester, and consists of 360 equiangular
projections over a full rotation.\
 The dataset has mechanical offsets and high noise,
with scan parameters set up for head CBCT i.e. low intensity X-rays (exact parameters unknown).

In this test, the sample has been reconstructed using FDK, CGLS, OS-SART and OS-ASD-POCS with 360 projections 
and 90 projections.

![im1](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/randofull.png)

RANDO head using 360 equiangular projections, with 4 different algorithms. From left to right, top to bottom: FDK, CGLS, OS-SART and OS-ASD-POCS. The displaying window is [0-0.05]

![im1](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/rando90.png)

RANDO head using 90 equiangular projections, with 4 different algorithms. From left to right, top to bottom: FDK, CGLS, OS-SART and OS-ASD-POCS. The displaying window is [0-0.05]

Special thanks to The Christie Hospital (Manchester) for sharing this dataset.

## Cryo soft X-rays

**Source of image** Ander Biguri's PhD thesis, page 101.

**Description:** This dataset contains a section of an image containing a few Trypanosoma, a unicellular
parasitic protozoa that cause different illnesses, such as the sleeping sickness.\
 Using soft X-ray on cryogenically frozen samples in synchrotron light sources, nanometre scale images can be achieved.\
 However due to the nature of the X-rays, the scans are very noisy and have very limited scan angle (less than 120 degrees).\
 This dataset has been acquired at the b24 beamline of the Diamond Light Source.
 
 In this experiment, the sample has been reconstructed using FDK, OS-SART and OS-ASD-POCS.
 
 ![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TV.png)\
 ![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TVz1.png)\
 ![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TVz2.png)

Columns: FBP, OS-SART (20 iterations) and OS-ASD-POCS (20 iterations,
20 TV iterations each). The red squares in the first row shows the location of the
zoomed-in areas from the second and third row.

Special thanks to the Diamond Light Source, specially Michele Darrow and Imanol Luengo for sharing this dataset.

## Conventional micro-CT

**Source of image:** Ander Biguri's PhD thesis, page 96.

**Description:** [SophiaBeads Dataset](https://zenodo.org/record/16474) are acquired specifically for testing and comparing reconstruction methods for X-ray computed tomography.\
The sample is a plastic tube with a diameter of 25 mm, filled with uniform Soda-Lime Glass (SiO2-Na2O) beads of diameters 2.5 mm (with standard deviation 0.1 mm). [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16474.svg)](https://doi.org/10.5281/zenodo.16474)\
The dataset is reconstructed with FDK and CGLS and the image profile shown over an arbitrary path.

![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/sophiaFDKCGLS.png)

SophiaBeads dataset with FDK and CGLS (15 iterations), using 256
projections. The red line shows the profile evaluated. The display window is [0-0.15].

![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/spophiaprofile.png)

Image profile on the SophiaBeads dataset with FDK and CGLS (15
iterations), using 256 projections

## Conventional micro-CT

**Source of image:** Arbitrarily large iterative tomographic reconstruction on multiple GPUs using the TIGRE toolbox [https://arxiv.org/abs/1905.03748](https://arxiv.org/abs/1905.03748) Fig 10.

**Description** Reconstruction of a panel shift micro-CT scan of a coffee bean. The reconstructed images are 3340 x 3340 x 900 voxels, with 2134 projections of size 3780 x 900. This reconstruction is to showcase possibility of very large reconstruction using TIGRE. Reconstruction using FDK and GCLS (a very memory heavy algorithm) are shown. 

**Images are too large for GitHub, click on the image description links to dowloand them**

![im](https://raw.githubusercontent.com/CERN/TIGRE/master/Frontispiece/coffee.png)

(Right) Coffee bean with [FDK](https://github.com/CERN/TIGRE/blob/master/Frontispiece/fdk.png?raw=true), (left) [CGLS](https://github.com/CERN/TIGRE/blob/master/Frontispiece/cgls.png?raw=true) 15 iterations, 2134 projections. *Image shown in very compressed format, please dowload originals for accurate representation.*






