TIGRE: Tomographic Iterative GPU-based Reconstruction Toolbox
======

TIGRE is a MATLAB/CUDA toolbox for fast and accurate 3D tomographic 
reconstruction. The aim of TIGRE is providing a wide range of easy-to-use 
iterative algorithms to the tomographic research community. 
We would like to build an stronger bridge between the algorithm researchers and
imaging researchers, by encouraging and supporting contributions frmo both sides into
TIGRE. TIGRE is free to use and distribute, use it, modify it, break it, share it; 
the only requirements is proper referencing to the authors.

Currently it contains Cone Beam CT geometries, and a beta of parallel geometires. 


**TIGRE features**:

  - Fast (state of the art) projection using interpolation or Siddon (Jacob) ray-tracing

  - Fast (state of the art) backprojection, with FDK weigth or matched weigths

  - A wide range of algorithms with multiple configurations for each in each
      - FDK                    
      - SART                    
      - OS-SART                
      - SIRT                   
      - CGLS
      - ADS-POCS               
      - OSC-POCS              
      - B-ADS-POCS-&#946;       
      - SART-TV

  - TV denosing for 3D images

  - A variety of plotting function

  - Quality measures


## How to install TIGRE

(Cauntion, this has only been tested in win64 machines, please [report any 
issue][2] if it doesnt work in other arch/SO)
 
   - Download TIGRE from the downloads page
   
   - Install  CUDA Toolkit (the later, the better)
     Download [here][1]
   
   - If you are working in a win64 machine, with CUDA 7.5, stop here. The
     toolbox should be ready to use.

If it doesnt work, or you are not in win64 CUDA 7.5
   
   - Make sure a compiler is set up in your MATLAB. run `mex -setup`. If a 
     compiler is not set up, make sure you download one (some are free)
     and run mex -setup again.
   
   - Make sure the xml file for compiling is proerly set up. E.g. in Linux 64 bit machines
     mex_CUDA_glnxa64.xml should be present, and inside it the proper links to CUDA have to be set up.

   - Run Compile.m

   - TIGRE should be ready to use!

## Getting started


The first thing you need to do is run `InitTIGRE`, to initialize all the 
folders and subfolders if the toolbox.

Currently the documentation is included in each of the functions. You can access it 
by typing `help function_name` or by selecting the function in the editor and pressing F1.
Additionally, the demos should include all the necessary documentation and examples of use.

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

**Q: I get "the launch timed out and was terminated" error when I run big images
in my computer**

*A: This happens because your GPU takes too long (according to the OS) to finish
running the code. Don't worry, too long means about 100ms or so. However, you need
to make sure to change the OS's GPU watchdog time. 
If you are working on a TESLA, setting the TESLA to TCC mode will fix the problem.*

**Q: After running something I got an error in Ax or Atb and now nothing works**

*A: Unfortunatedly when CUDA has an error, it hungs. You need to restart MATLAB to fix
this. Hopefully we can find a solution for this without the need of restarting MATLAB*

## Licensing

TIGRE toolbox is released under the BSD License, meaning you can use and modify 
the software freely in any case, but you **must** give attribution to the original authors.
For more information, read the license file or the [BSD License Definition][5] or [the license file][6]

## Contact

*Before contacting* consider that you migth be able to let us know your concerns by
[raising an issue][2] and [labeling accordingly][4]. If you still need to contact us:

tigre.toolbox@gmail.com

We will make an effort to answer as soon as we can.

## Referencing TIGRE

If you use TIGRE in any publications, please reference the following paper:

**TIGRE: A MATLAB-GPU toolbox for CBCT image reconstruction**
*Ander Biguri, Manjit Dosanjh, Steven Hancock, and Manuchehr Soleimani*
*To be published*, 2016

Also, if you use any algorithm, you should reference the corresponding creator
of the algorithms. If you dont know the article, use `citeme('NameOfFunction')`
and the rigth reference will appear.


---

<sup>1</sup> You can learn, it is very usefull!


[1]: https://developer.nvidia.com/cuda-downloads
[2]: https://github.com/CERN/TIGRE/issues
[3]: https://stackoverflow.com/help/mcve
[4]: https://help.github.com/articles/applying-labels-to-issues-and-pull-requests/
[5]: http://www.linfo.org/bsdlicense.html
[6]: https://github.com/CERN/TIGRE/license.txt
