FAQ
======

If you have any issues that is not addressed in the FAQ, or you found a bug in
the code, please report it on [the issues page][2].

If it is a specific problem in a specific scenario, please provide a [*Mnimal Complete Verifiable Example*][3] of the problem.
Make it as short as possible. Ensure that your example is self-complete, no other code is needed. Ensure that someone else can execute it and 
reproduce the error you are having. 


***

**Q: I can not compile TIGRE!**

*A: Have a look at the installation instructions of [the MATLAB version](Frontispiece/MATLAB_installation) or [the Python version](Frontispiece/Python_installation).\
If after reading that you still have problems, please feel free to contact us.*

**Q: TIGRE succesfully compiles but I get `texture memory fail` when I run the code**

*A: It is likely you are running on a GPU that has compute capability older than 3.0, which is not supported.\
 However, you may be able to still run TIGRE, by compiling with CUDA 8.0 and changing line 38 of `mex_CUDA_win.xml` to:*
 ```
 COMPFLAGS="-gencode=arch=compute_21,code=sm_21 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_50,code=&#92;&quot;sm_50,compute_50&#92;&quot; -gencode arch=compute_52,code=sm_52 --ptxas-options=-v --compiler-options=/c,/GR,/W3,/EHs,/nologo,/MD"
 ```

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

**Q: Is the pyhton versiona vailable for Windows?**

*A: Not yet. We do want it running but we have problems with compilation. Feel free to contribute, we'd like help with this problem*

**Q: Can I use your code for XXXXXX**

*A: Yes. The BSD license is a fully permisive one. Your onyl requirement is to reference the usage of TIGRE, or that you based your work on TIGRE*

**Q: I want to collaborate, how can I do it?**

*A: We are very happy to hear that! Please do [contact us](mailto:ander.biguri@gmail.com), or use the standard Github mechanisms to add/fix code*

[2]: https://github.com/CERN/TIGRE/issues
[3]: https://stackoverflow.com/help/mcve