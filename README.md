<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-13-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->


<p align="center">
  <a href="https://cern.github.io/TIGRE/">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="docs/_static/logo/tigre-logo-dark.svg">
      <img alt="TIGRE" src="docs/_static/logo/tigre-logo.svg" width="330">
    </picture>
  </a>
</p>

TIGRE: Tomographic Iterative GPU-based Reconstruction Toolbox
======

TIGRE is an open-source toolbox for fast and accurate 3D tomographic 
reconstruction for any geometry.  Its focus is on iterative algorithms 
for improved image quality that have all been optimized to run on GPUs 
(including multi-GPUs) for improved speed. It combines the higher level 
abstraction of MATLAB or Python with the performance of CUDA at a lower level in order to make 
it both fast and easy to use.

**Website and documentation: https://cern.github.io/TIGRE/**

Read the most up to date article about tigre [here](https://iopscience.iop.org/article/10.1088/2631-8695/adbb3a).

**[Website](https://cern.github.io/TIGRE/)** ·
**[Install](https://cern.github.io/TIGRE/install/index.html)** ·
**[User guide](https://cern.github.io/TIGRE/user-guide/index.html)** ·
**[Demos](https://cern.github.io/TIGRE/demos/index.html)** ·
**[API reference](https://cern.github.io/TIGRE/api/index.html)** ·
**[Gallery](https://cern.github.io/TIGRE/gallery/index.html)** ·
**[Cite TIGRE](https://cern.github.io/TIGRE/citing.html)** ·
**[Release notes](https://cern.github.io/TIGRE/changelog.html)**

TIGRE is free to download and distribute: use it, modify it, add to it,
share it. Our aim is to provide a wide range of easy-to-use algorithms
for the tomographic community "off the shelf", and to build a bridge between
algorithm developers and imaging researchers and clinicians. Requests and
ideas are welcome in the [discussions](https://github.com/CERN/TIGRE/discussions).

## Features

- **MATLAB** and **Python**, with the same CUDA code, algorithms and demos.
- Fast projection and backprojection on **GPUs**, including **multiple GPUs**
  and images larger than GPU memory.
- **Any geometry**: cone beam, parallel beam, fan beam, helical, tomosynthesis,
  C-arm, curved detectors. Geometric parameters are defined per projection.
- **Algorithms**: FDK/FBP; the SART family (SART, OS-SART, SIRT); total-variation
  methods (ASD-POCS and its variants, SART-TV, FISTA); Krylov methods (CGLS,
  LSQR, hybrid LSQR, LSMR, IRN-TV-CGLS, hybrid fLSQR-TV, AB/BA-GMRES); MLEM;
  prior-image PICCS.
- **Real data**: loaders for Nikon, Bruker, YXLON, Diondo and Varian scanners,
  DXchange (synchrotron) and DICOM.
- PyTorch operators with gradients, plotting, image quality metrics and TV
  denoising.

See the [user guide](https://cern.github.io/TIGRE/user-guide/index.html) for
how to choose an algorithm, and the
[Gallery](https://cern.github.io/TIGRE/gallery/index.html) for reconstructions
of real scans.

## Installation

> [!IMPORTANT]
> TIGRE is **not** on PyPI. `pip install tigre` installs an unrelated package with the same name.

- **Python on Linux:** prebuilt wheels, no compiler needed. Copy the command for
  your Python version from the [Python install page](https://cern.github.io/TIGRE/install/python.html).
- **Python on Windows, or from source:** `pip install git+https://github.com/CERN/TIGRE.git`
  (needs the CUDA toolkit and a C++ compiler; see the
  [Python install page](https://cern.github.io/TIGRE/install/python.html)).
- **MATLAB:** run `Compile` in the `MATLAB` folder; see the
  [MATLAB install page](https://cern.github.io/TIGRE/install/matlab.html).

Then start with the [demos](https://cern.github.io/TIGRE/demos/index.html)
(also in [`Python/demos`](Python/demos) and [`MATLAB/Demos`](MATLAB/Demos)).
If something goes wrong, see
[Troubleshooting and FAQ](https://cern.github.io/TIGRE/user-guide/troubleshooting.html).

## Contact

Contact the authors directly at:

[tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com) or [ander.biguri@gmail.com](mailto:ander.biguri@gmail.com)

for any questions/comments or if you want to be added to the mailing list or the Slack team.

The Slack team is a good place for chatting about development and questions about TIGRE. Please send an email to the authors and you will receive an invitation.

## Licensing

The creation of TIGRE was supported by the University of Bath and CERN. It is released under the BSD License, meaning you can use and modify the software freely.  However, you **must** cite the original authors.
For more information read [the license file][1] or the [BSD License Definition][2].

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

You can find the most up to date article in the TIGREv3 release paper, with most novel features well documented.

**TIGRE v3: Efficient and easy to use iterative computed tomographic reconstruction toolbox for real datasets**
*Ander Biguri, Tomoyuki Sadakane, Reuben Lindroos, Yi Liu, Malena Sabaté Landman, Yi Du, Manasavee Lohvithee, Stefanie Kaser, 
Sepideh Hatamikia, Robert Bryll, Emilien Valat, Sarinrat Wonglee, Thomas Blumensath and Carola-Bibiane Schönlieb*
**Engineering Research Express**

[Read the article (open access)][6]

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
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/malena-sabate"><img src="https://avatars.githubusercontent.com/u/83228154?v=4?s=100" width="100px;" alt="malena-sabate"/><br /><sub><b>malena-sabate</b></sub></a><br /><a href="#research-malena-sabate" title="Research">🔬</a> <a href="https://github.com/CERN/TIGRE/commits?author=malena-sabate" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://emilienvalat.net/"><img src="https://avatars.githubusercontent.com/u/46785587?v=4?s=100" width="100px;" alt="Emilien Valat"/><br /><sub><b>Emilien Valat</b></sub></a><br /><a href="#research-Emvlt" title="Research">🔬</a> <a href="https://github.com/CERN/TIGRE/commits?author=Emvlt" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/zezisme"><img src="https://avatars.githubusercontent.com/u/40331996?v=4?s=100" width="100px;" alt="zez666"/><br /><sub><b>zez666</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/issues?q=author%3Azezisme" title="Bug reports">🐛</a> <a href="https://github.com/CERN/TIGRE/commits?author=zezisme" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Pawelip"><img src="https://avatars.githubusercontent.com/u/84121133?v=4?s=100" width="100px;" alt="Pawel_Lipowicz"/><br /><sub><b>Pawel_Lipowicz</b></sub></a><br /><a href="https://github.com/CERN/TIGRE/commits?author=Pawelip" title="Code">💻</a> <a href="https://github.com/CERN/TIGRE/commits?author=Pawelip" title="Tests">⚠️</a> <a href="#userTesting-Pawelip" title="User Testing">📓</a></td>
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
[6]: https://iopscience.iop.org/article/10.1088/2631-8695/adbb3a
