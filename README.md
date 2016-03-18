TIGRE: Tomographic Iterative GPU-based Reconstruction Toolbox
======

TIGRE is a MATLAB/CUDA toolbox for fast and accurate 3D tomographic 
reconstruction. Currently it contains Cone Beam CT geometries. The
objective of TIGRE is to provide tomograpy users a wide variety of image 
reconstruction algorithims in an easy to use and fast form, and image 
reconstruction reearchers basic blocks for algorithm design and testing for
big datasets.

**TIGRE features**:

  - Fast (state of the art) projection using interpolation or Siddon ray-tracing

  - Fast (state of the art) backprojection, with FDK weigth or matched weigths

  - A wide range of algorithms with multiple parameters in each
      - FDK                    
      - SART                    
      - OS-SART                
      - SIRT                   
      - CGLS
      - ADS-POCS               
      - OSC-POCS              
      - B-ADS-POCS-<p>&beta<p>         
      - SART-TV

  - TV denosing for 3D images

  - A variety of plotting function

  - Quality measures


## How to install TIGRE

(Cauntion, this has only been tested in win64 machines, please [report any 
issue][2] if it doesnt work in other arch/SO)
 
   - Download TIGRE from the downloads page
   
   - Install  CUDA Toolit (the later, the better)
     Download [here][1]
   
   - If you are working in a win64 machine, with CUDA 7.5, stop here. The
     toolbox should be ready to use.

If it doesnt work, or you are not in win64 CUDA 7.5
   
   - Make sure a compiler is set up in your MATLAB. run `mex -setup`. If a 
     compiler is not set up, make sure you download one (some are free)
     and run mex -setup again.

   - Run Compile.m

   - TIGRE should be ready to use!

## Getting started


The first thing you need to do is run `InitTIGRE`, to initialize all the 
folders and subfolders if the toolbox.

We are currently working in creating some good documentation for TIGRE, but
meanwhile, you can go to the demos folder. There are several examples on everything
you need to run TIGRE.

## Issues

If you have any issues with compiling/running TIGRE, or you found a bug in
the code, please report it on [the issues page][2]<sup>*</sup>.

Please, be as clear as  possible with the description of the problem.
If it is a specific problem in an specific scenario, provide a *Minimum 
Complete Verifiable Example* of the problem (see [Stackoverlfowcom's definition][3]);

<sup>*</sup> **Protip** you can [label the Issues by topic][4]

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

## Licensing


## Contact

*Before contacting* consider that you migth be able to let us know your concerns by
[raising an issue][2] and [labeling accordingly][4]. If you still need to contact us:

tigre.toolbox@gmail.com

We will make an effort to answer as soon as we can.

## Referencing TIGRE

If you use TIGRE in any publications, please reference the following paper:



Also, if you use any algorithm, you should reference the corresponding creator
of the algorithms. If you dont know the article, use `citeme('NameOfFunction')`
and the rigth reference will appear.


---

<sup>1</sup> You can learn, it is very usefull!


[1]: https://developer.nvidia.com/cuda-downloads
[2]: https://github.com/AnderBiguri/TIGRE/issues
[3]: https://stackoverflow.com/help/mcve
[4]: https://help.github.com/articles/applying-labels-to-issues-and-pull-requests/