TIGRE: Tomographic Iterative GPU-based Reconstruction Toolbox
======

TIGRE is a PYTHON/CUDA toolbox for fast and accurate 3D tomographic 
reconstruction, created jointly by University of Bath's Engineering Tomography Lab and CERN. 
The aim of TIGRE is to provide a wide range of easy-to-use iterative algorithms to the tomographic research community. 
We would like to build an stronger bridge between algorithm researchers and
imaging researchers, by encouraging and supporting contributions from both sides into
TIGRE. TIGRE is free to use and distribute, use it, modify it, break it, share it; 
the only requirement is proper referencing to the authors.

Currently it contains Cone Beam CT geometries, and a beta of parallel geometries. 

If you wissh to be added to the mailing list of TIGRE, please send an email to tigre.toolbox@gmail.com. 
We will only contact you whenever a new version of TIGRE is released or with important news.

TIGRE is still being developed and we are still adding new features to it. If you have any request for your specific application, do not hesitate to contact us!

**TIGRE in media**:

[Article of Medical Physics Web on TIGRE](http://medicalphysicsweb.org/cws/article/research/66343)



Donwload: [https://github.com/CERN/TIGRE/releases][7]


**NEWS**

TIGRE-python version 0.1.0 released for use. Changes have been made to the architecture to make it more modular. 
As well as this, style follows the formal python format more closely; classes are either all upper case or
SemiUpperCase. Functions, such as the algorithms, are lowercase and are extensions of the umbrella class
IterativeReconAlg. Demos have been updated for this format and can be launched from the ipython console. More detailed
instructions will follow as progress continues towards making a more stable and modular TIGRE-python.


**TIGRE features**:

  - Fast, state of the art projection using interpolation or Siddon (Jacob) ray-tracing

  - Fast, state of the art backprojection, with FDK weight or matched weight 

  - A wide range of algorithms with multiple configurations for each 
      - FDK
      - Gradien descend family
        - SART                    
        - OS-SART                
        - SIRT

  - A variety of plotting functions

  - Quality measures

## How to install TIGRE-python

(Tested on Linux 64 machines, please [report any 
issue][2] if it doesnt work in other arch/OS)

   - Install  CUDA Toolkit (the later, the better)
     Download [here][1]
   - Clone the repo by:  
      > git clone https://github.com/CERN/TIGRE.git
   - Change directory 
      > cd TIGRE/Python/ 
   - and run the setup file:  
      > python setup.py install
   - python may need permission to install to certain directories, 
   in which case do:
      > sudo --preserve-env python setup.py install

If it doesn't work

   - Try setting up a py27 environment with: 

       > conda create -n py27 python=2.7 anaconda

   - Now activate the environment from the CLI with 

       > source activate py27
 
   - try the above steps again (making sure the py27 environment is activated)
 
   - TIGRE should be ready to use!

**Getting started with TIGRE-python**

  - After installing the software following the instructions below, we strongly reccomend downloading the tigre demonstration file from the   repo and looking at the example code there. The files are in jupyter notebook form and should be self contained (as long as the software is correctly installed). 

## Issues

If you have any issues with compiling/running TIGRE, or you found a bug in
the code, please report it on [the issues page][2]<sup>&#8727;</sup>.

Please, be as clear as  possible with the description of the problem.
If it is a specific problem in an specific scenario, provide a *Minimum 
Complete Verifiable Example* of the problem (see [Stackoverflow.com's definition][3]);

<sup>&#8727;</sup> **Protip** you can [label the Issues by topic][4]

## I want to contribute! 

Wow! We are glad you do! Please, if you want to contribute new methods, 
algorithms, pre- or post- processing techniques, bug-fixes, documentation, or
*anything* you thing it can help the community, please, do! We encourage
this behaviour!

*But how?* 

Easy! you can download the git repo, and make a pull request with your 
contribution. We will review it and add it if its suited for TIGRE. 

If you don't know how git
works<sup>1</sup> you can also send an email us to tigre.toolbox@gmail.com 
with your contribution, and an explanation of what it does and how. We will
review and add the contribution if suited.

If your contribution can be linked to a published, peer-reviewed article or
an arXiv entry, please let us know so we can make sure the `citeme` function
includes your contributions.

## FAQ

**Q: I get a fair amount of warnings when I compile the code, what is happening?**

*A: Do not worry about the warnings. We are perfectly aware of them and know that they will have no 
effect whatsoever in the correct execution (both in computational time and accuracy) of the code.*

## Licensing

TIGRE toolbox is released under the BSD License, meaning you can use and modify 
the software freely in any case, but you **must** give attribution to the original authors.
For more information, read the license file or the [BSD License Definition][5] or [the license file][6]

## Contact

*Before contacting* consider that you might be able to let us know any problem TIGRE may have by
[raising an issue][2] and [labeling accordingly][4]. 

If you want to contact us for other reasons than an issue with the tool, please send us an email to

tigre.toolbox@gmail.com

or contact the author directly in

MATLAB: 
a.biguri@bath.ac.uk

Python:
reuben.jonathan.lindroos@cern.ch 

We will make an effort to answer as soon as we can.

## Referencing TIGRE

If you use TIGRE in any publications, please reference the following paper:

**TIGRE: A MATLAB-GPU toolbox for CBCT image reconstruction**
*Ander Biguri, Manjit Dosanjh, Steven Hancock, and Manuchehr Soleimani*
**Biomedical Physics & Engineering Express, Volume 2, Number 5**
[Read the article (Open Access)][8]

Also, if you use any algorithm, you should reference the corresponding creator
of the algorithms. If you don't know the article, use `citeme('NameOfFunction')`
and the right reference will appear.


---

<sup>1</sup> You can learn, it is very useful!


[1]: https://developer.nvidia.com/cuda-downloads
[2]: https://github.com/CERN/TIGRE/issues
[3]: https://stackoverflow.com/help/mcve
[4]: https://help.github.com/articles/applying-labels-to-issues-and-pull-requests/
[5]: http://www.linfo.org/bsdlicense.html
[6]: https://github.com/CERN/TIGRE/license.txt
[7]: https://github.com/CERN/TIGRE/releases
[8]: http://iopscience.iop.org/article/10.1088/2057-1976/2/5/055010
