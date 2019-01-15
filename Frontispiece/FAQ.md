FAQ
======

If you have any issues that is not addressed in the FAQ, or you found a bug in
the code, please report it on [the issues page][2]<sup>&#8727;</sup>.

If it is a specific problem in a specific scenario, please provide a [*Minimum 
Complete Verifiable Example*][3] of the problem.




**Q: I get a fair amount of warnings when I compile the code, what is happening?**

*A: Do not worry about the warnings. We are perfectly aware of them and know that they will have no 
effect whatsoever in the correct execution (both in computational time and accuracy) of the code.*

**Q: I get "the launch timed out and was terminated" error when I run big images
in my computer**

*A: This happens because your GPU takes too long (according to the OS) to finish
running the code. Don't worry, too long means about 100ms or so. However, you need
to make sure to change the OS's GPU watchdog time. 
If you are working on a TESLA, setting the TESLA to TCC mode will fix the problem.*

**Q: After running something I got an error in Ax or Atb and now nothing works**

*A: Unfortunately when CUDA has an error, it hungs. You need to restart MATLAB to fix
this. Hopefully we can find a solution for this without the need of restarting MATLAB*

**Q: Does it work in MATLAB XXXX with Compiler XXXX in OS XXXX**

*A: In general, it should, as long as you are following both CUDAs and MATLABs s upportedcompilers.
The truth is that there is few compilers that fit the eeligibility criteria for both. MATLAB version and OS 
should not be a problem.*



[2]: https://github.com/CERN/TIGRE/issues
[3]: https://stackoverflow.com/help/mcve