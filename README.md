[![Documentation Status](https://readthedocs.org/projects/tigre/badge/?version=latest)](https://tigre.readthedocs.io/en/latest/?badge=latest)
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-9-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->


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
application, do not hesitate to [contact us](#contact) or open a  [discussion thread](https://github.com/CERN/TIGRE/discussions)!

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
       
        - Gradient-based algorithms (SART, OS-SART, SIRT, ASD-POCS, OS-ASD-POCS, B-ASD-POCS-β, PCSD, AwPCSD, Aw-ASD-POCS) with multiple tuning parameters (Nesterov acceleration, initialization, parameter reduction, ...)
       
        - Krylov subspace algorithms (CGLS, LSQR, hybrid LSQR, LSMR, IRN-TV-CGLS, hybrid-fLSQR-TV, AB/BA-GMRES)
       
        - Statistical reconstruction (MLEM)
       
        - Variational methods (FISTA, SART-TV) 
       
- TV denoising for 3D images.
       
- Basic image loading functionality.
       
- A variety of plotting functions.
       
- Image quality metrics.

- Nikon and Varian and Phillips (DICOM) scanner data loaders. 

## Installation

MATLAB and Python builds are both fully supported.

- [Installation instructions and requirements for MATLAB](Frontispiece/MATLAB_installation.md).

- [Installation instructions and requirements for Python](Frontispiece/python_installation.md). 

**Advanced, not required to run TIGRE**, will change the source code. Only do if performance is critical.

- [Tune TIGRE for machine. Tricks to slightly speed up the code](Frontispiece/Tune_TIGRE.md)


## FAQ

For answers to frequently asked questions [click here](Frontispiece/FAQ.md).

If you have new question not answered in the FAQ, please [contact us](#contact), join the [Slack group](#contact) or open a  [discussion thread](https://github.com/CERN/TIGRE/discussions).

## Gallery

To see a gallery of images of different CT modalities reconstructed using TIGRE [click here](Frontispiece/Gallery.md).

<img src="https://raw.githubusercontent.com/AnderBiguri/PhDThesis/master/Applications/randofull.png" height="400">


## Further Reading

If you want more information on TIGRE and its algorithms, [click here](Frontispiece/Further_reading.md).


## Contact

Contact the authors directly at:

[tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com) or [ander.biguri@gmail.com](mailto:ander.biguri@gmail.com)

for any questions/comments or if you want to be added to the mailing list or the Slack team.

The Slack team is a good place for chatting about development and questions about TIGRE. Please send an email to the authors and you will receive an invitation.

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

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/AnderBiguri"><img src="https://avatars.githubusercontent.com/u/11854388?v=4?s=100" width="100px;" alt="Biguri"/><br /><sub><b>Biguri</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=AnderBiguri" title="Code">💻</a> <a href="#example-AnderBiguri" title="Examples">💡</a> <a href="#ideas-AnderBiguri" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-AnderBiguri" title="Maintenance">🚧</a> <a href="#research-AnderBiguri" title="Research">🔬</a> <a href="https://github.com/CERN/TIGRE/pulls?q=is%3Apr+reviewed-by%3AAnderBiguri" title="Reviewed Pull Requests">👀</a> <a href="#tutorial-AnderBiguri" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/yliu88au"><img src="https://avatars.githubusercontent.com/u/75292881?v=4?s=100" width="100px;" alt="yliu88au"/><br /><sub><b>yliu88au</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=yliu88au" title="Code">💻</a> <a href="https://github.com/CERN/TIGRE/issues?q=author%3Ayliu88au" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/reubenlindroos"><img src="https://avatars.githubusercontent.com/u/25688713?v=4?s=100" width="100px;" alt="Reuben Lindroos"/><br /><sub><b>Reuben Lindroos</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=reubenlindroos" title="Code">💻</a> <a href="https://github.com/CERN/TIGRE/issues?q=author%3Areubenlindroos" title="Bug reports">🐛</a> <a href="#design-reubenlindroos" title="Design">🎨</a> <a href="#ideas-reubenlindroos" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/genusn"><img src="https://avatars.githubusercontent.com/u/25704789?v=4?s=100" width="100px;" alt="genusn"/><br /><sub><b>genusn</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=genusn" title="Code">💻</a> <a href="https://github.com/CERN/TIGRE/issues?q=author%3Agenusn" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/yidu-bjcancer"><img src="https://avatars.githubusercontent.com/u/7495679?v=4?s=100" width="100px;" alt="Yi DU"/><br /><sub><b>Yi DU</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=yidu-bjcancer" title="Code">💻</a> <a href="https://github.com/CERN/TIGRE/issues?q=author%3Ayidu-bjcancer" title="Bug reports">🐛</a> <a href="#research-yidu-bjcancer" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/tsadakane"><img src="https://avatars.githubusercontent.com/u/40597344?v=4?s=100" width="100px;" alt="tsadakane"/><br /><sub><b>tsadakane</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=tsadakane" title="Code">💻</a> <a href="#design-tsadakane" title="Design">🎨</a> <a href="#ideas-tsadakane" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/CERN/TIGRE/issues?q=author%3Atsadakane" title="Bug reports">🐛</a> <a href="#tutorial-tsadakane" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://omg.lol/sco1"><img src="https://avatars.githubusercontent.com/u/5323929?v=4?s=100" width="100px;" alt="S. Co1"/><br /><sub><b>S. Co1</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=sco1" title="Code">💻</a> <a href="#design-sco1" title="Design">🎨</a> <a href="#tool-sco1" title="Tools">🔧</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Daveelvt"><img src="https://avatars.githubusercontent.com/u/16086944?v=4?s=100" width="100px;" alt="Daveelvt"/><br /><sub><b>Daveelvt</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=Daveelvt" title="Code">💻</a> <a href="#research-Daveelvt" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/phernst"><img src="https://avatars.githubusercontent.com/u/9623894?v=4?s=100" width="100px;" alt="phernst"/><br /><sub><b>phernst</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/issues?q=author%3Aphernst" title="Bug reports">🐛</a> <a href="https://github.com/CERN/TIGRE/commits?author=phernst" title="Code">💻</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->


[1]: LICENSE.txt
[2]: http://www.linfo.org/bsdlicense.html
[3]: http://iopscience.iop.org/article/10.1088/2057-1976/2/5/055010
[4]: https://doi.org/10.1016/j.jpdc.2020.07.004
[5]: https://arxiv.org/abs/1905.03748
