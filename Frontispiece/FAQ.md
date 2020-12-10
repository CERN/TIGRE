FAQ
======

If you have any issues that is not addressed in the FAQ, or you found a bug in
the code, please report it on [the issues page][2].

If it is a specific problem in a specific scenario, please provide a [*Minimal Complete Verifiable Example*][3] of the problem.
Make it as short as possible. Ensure that your example is self-complete, no other code is needed. Ensure that someone else can execute it and 
reproduce the error you are having. 

If you have more of a comment or a suggestion, you can also open a thread in the  [discussions](https://github.com/CERN/TIGRE/discussions) page. 

***

**Q: I can not compile TIGRE!**

*A: Have a look at the installation instructions of [the MATLAB version](Frontispiece/MATLAB_installation) or [the Python version](Frontispiece/Python_installation).\
If after reading that you still have problems, please feel free to contact us.*

**Q: I get compilation error nvcc fatal   : Unsupported gpu architecture 'compute_XX'**

Your particular CUDA version does not support all compute capablities of all GPUs. TIGRE by default compiles assuming
CUDA 11.0, meaning we compile up to compute capability 8.6 (`-gencode=arch=compute_86,code=sm_86`). 
If you are using an older CUDA version, maybe you can not compile to some architectures. Try removing the higher numbered arguments in the following line in the XML,
or a similar line in the `setup.py` file for python. 

 ```
COMPFLAGS="-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_86,code=sm_86 --default-stream per-thread  --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"
 ```

**Q: TIGRE succesfully compiles but I get `texture memory fail` when I run the code**

*A: It is likely you are running on a GPU that has compute capability older than 3.5, which is not supported.\

**Q: I get a fair amount of warnings when I compile the code, what is happening?**

*A: Do not worry about the warnings.\
We are perfectly aware of them and know that they will have no 
effect whatsoever in the correct execution (both in computational time and accuracy) of the code.*

**Q: I get "the launch timed out and was terminated" error when I run big images
in my computer**

*A: This happens because your GPU takes too long (according to the OS) to finish
running the code. Don't worry, too long means about 100ms or so. However, you need
to make sure to change the OS's GPU watchdog time. 
If you are working on a TESLA, setting the TESLA to TCC mode will fix the problem.*

**Q: I get "unespecified kernel launch failure error when I run big images
in my computer**

*A: This also happens because your GPU takes too long (according to the OS) to finish
running the code. Don't worry, too long means about 100ms or so. However, you need
to make sure to change the OS's GPU watchdog time. \
Go to `NVIDIA Nsight->Nsigth monitor options->WDDM TDR Enabled->False`\
Alternatively set `NVIDIA Nsight->Nsigth monitor options->WDDM TDR Delay` to a higher number*

**Q: After running something I got an error in Ax or Atb and now nothing works**

*A: Unfortunately when CUDA has an error, it hungs. You need to restart MATLAB to fix
this. Hopefully we can find a solution for this without the need of restarting MATLAB*

**Q: Does it work in Windows with MinGW**

*A: No. CUDA only allows windows compilation using Visual Studio. The free versions are good enough, however.*

**Q: Is the pyhton version available for Windows?**

*A: Yes, but there are particular versions of python 3 that had caused problems, so please report back to us if you had issues. *

**Q: Can I use your code for XXXXXX**

*A: Yes. The BSD license is a fully permisive one. Your only requirement is to reference the usage of TIGRE, or that you based your work on TIGRE*

**Q: I want to collaborate, how can I do it?**

*A: We are very happy to hear that! Please do [contact us](mailto:ander.biguri@gmail.com), or use the standard Github mechanisms to add/fix code*

[2]: https://github.com/CERN/TIGRE/issues
[3]: https://stackoverflow.com/help/mcve