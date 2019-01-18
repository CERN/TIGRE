Gallery
=======

This page showcases some images of applications of TIGRE and directs to papers using TIGRE in a variety of CT cases.


## Few images of TIGRE algorithms in different CT applications. 

### RANDO Head

**Source of image:** Ander Biguri's PhD thesis, page 93.

**Description:** This dataset is a RANDO head phantom that represents a head of an adult. The
dataset was acquired in The Christie Hospital in Manchester, and consists of 360 equiangular
projections over a full rotation.\
 The dataset has mechanical offsets and high noise,
with scan parameters set up for head CBCT i.e. low intensity X-rays (exact parameter sunknown).

In this test, the sample has been reconstructed using FDK, CGLS, OS-SART and OS-ASD-POCS with 360 projections 
and 90 projections.

![im1](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/randofull.png)

RANDO head using 360 equiangular projections, with 4 different algorithms. From left to right, top to bottom: FDK, CGLS, OS-SART and OS-ASD-POCS. The displaying window is [0-0.05]

![im1](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/rando90.png)

RANDO head using 90 equiangular projections, with 4 different algorithms. From left to right, top to bottom: FDK, CGLS, OS-SART and OS-ASD-POCS. The displaying window is [0-0.05]

Special thanks to The Christie Hospital (Manchester) for sharing this dataset.

### Cryo soft X-rays

**Source of image** Ander Biguri's PhD thesis, page 101.

**Description:** This dataset contains a section of an image containing a few Trypanosoma, a unicellular
parasitic protozoa that cause dfferent illnesses, such as the sleeping sickness.\
 Using soft X-ray on cryogenically frozen samples in syncrotron light sources, nanometer scale images can be achieved.\
 However due to the nature of the X-rays, the scans are very noisy and have very limited scan angle (less than 120 degrees).\
 This dataset has been adquired at the b24 beamline of the Diamond Light Source.
 
 In this experiment, the sample has been reconstructed using FDK, OS-SART and OS-ASD-POCS.
 
 ![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TV.png)\
 ![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TVz1.png)\
 ![im](https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/FBP_OSSART_TVz2.png)

Columns: FBP, OS-SART (20 iterations) and OS-ASD-POCS (20 iterations,
20 TV iterations each). The red squares in the first row show sthe location of the
zoomed-in areas from the second and third row.
