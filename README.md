[![Documentation Status](https://readthedocs.org/projects/tigre/badge/?version=latest)](https://tigre.readthedocs.io/en/latest/?badge=latest)


TIGRE: Tomographic Iterative GPU-based Reconstruction Toolbox
======

TIGRE is an open-source toolbox for fast and accurate 3D tomographic 
reconstruction for any geometry.  Its focus is on iterative algorithms 
for improved image quality that have all been optimized to run on GPUs 
(including multi-GPUs) for improved speed. It combines the higher level 
abstraction of MATLAB or Python with the performance of CUDA at a lower level in order to make 
it both fast and easy to use.

TIGRE is free to download and distribute: use it, modify it, add to it, 
share it. Our aim is to provide a wide range of easy-to-use algorithms 
for the tomographic community "off the shelf".  We would like to build a 
stronger bridge between algorithm developers and imaging 
researchers/clinicians by encouraging and supporting contributions from 
both sides into TIGRE.

TIGRE remains under development as we are still adding new features 
(e.g., motion compensation).  If you have any request for a specific 
application, do not hesitate to [contact us](#contact)!

 - [TIGRE features](#features)
 
 - [Installation instructions](#installation)
 
 - [FAQ](#faq)
  
 - [Further reading](#further-reading)
 
 - [Contact](#contact) 
 
 - [Licensing](#licensing)


## TIGRE features

TIGRE is a GPU-based CT reconstruction software repository that contains a wide variety of iterative algorithms.

- **MATLAB** and **Python** libraries for high-performance x-ray absorption tomographic reconstruction.

- State-of-the-art implementations of projection and backprojection operations on **GPUs** (including **multi-GPUs**), with a simple interface using higher level languages to facilitate the development of new methods.

- **Flexible CT geometry:** Cone Beam, Parallel Beam, Digital Tomosynthesis, C-arm CT, and any other geometry.  Geometric parameters are defined per projection, not per scan.

- A wide range of reconstruction algorithms for CT.

    - Filtered backprojection (FBP,FDK) and variations (different filters, Parker weights, ...)
   
    - **Iterative algorithms**
       
        - Gradient-based algorithms (SART, OS-SART, SIRT) with multiple tuning parameters (Nesterov acceleration, initialization, parameter reduction, ...)
       
        - Krylov subspace algorithms (CGLS)
       
        - Statistical reconstruction (MLEM)
       
        - Total variation regularization based algorithms: proximal-based (FISTA, SART-TV) and POCS-based (ASD-POCS, OS-ASD-POCS, B-ASD-POCS-β, PCSD, AwPCSD, Aw-ASD-POCS)
       
- TV denoising for 3D images.
       
- Basic image loading functionality.
       
- A variety of plotting functions.
       
- Image quality metrics.
    

## Installation

MATLAB and Python builds are both fully supported.

- [Installation instructions and requirements for MATLAB](Frontispiece/MATLAB_installation.md).

- [Installation instructions and requirements for Python](Frontispiece/python_installation.md). 

**Advanced, not required to run TIGRE**, will change the source code. Only do if performance is critical.

- [Tune TIGRE for machine. Tricks to slightly speed up the code](Frontispiece/Tune_TIGRE.md)


## FAQ

For answers to frequently asked questions [click here](Frontispiece/FAQ.md).


## Gallery

To see a gallery of images of different CT modalities reconstructed using TIGRE [click here](Frontispiece/Gallery.md).

<img src="https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/randofull.png" height="400">


## Further Reading

If you want more information on TIGRE and its algorithms, [click here](Frontispiece/Further_reading.md).


## Contact

Contact the authors directly at:

[tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com) or [ander.biguri@gmail.com](mailto:ander.biguri@gmail.com)


## Licensing

The creation of TIGRE was supported by the University of Bath and CERN. It is released under the BSD License, meaning you can use and modify the software freely.  However, you **must** cite the original authors.
For more information read [the licence file][1] or the [BSD License Definition][2].

If you use TIGRE, please reference the following papers:

**TIGRE: A MATLAB-GPU toolbox for CBCT image reconstruction**
*Ander Biguri, Manjit Dosanjh, Steven Hancock and Manuchehr Soleimani*
**Biomedical Physics & Engineering Express, Volume 2, Number 5**
[Read the article (open access)][3]

And especially if you use images bigger than 512<sup>3</sup> or multiple GPUs

**Arbitrarily large iterative tomographic reconstruction on multiple GPUs using the TIGRE toolbox**
*Ander Biguri, Reuben Lindroos, Robert Bryll, Hossein Towsyfyan, Hans Deyhle, Ibrahim El khalil Harrane, Richard
Boardman, Mark Mavrogordato, Manjit Dosanjh, Steven Hancock, Thomas Blumensath*
**Journal of Parallel and Distributed Computing**
[Read the article][4], 
[Preprint][5]

[1]: LICENSE.txt
[2]: http://www.linfo.org/bsdlicense.html
[3]: http://iopscience.iop.org/article/10.1088/2057-1976/2/5/055010
[4]: https://doi.org/10.1016/j.jpdc.2020.07.004
[5]: https://arxiv.org/abs/1905.03748
