Changelog
======

Please, add here all changes added since previous version of TIGRE (v2.1)

## Features

- You can now select which GPUs to run the code on, either by name, or GPUid. By @genusn
- Improved compilation, now there is no need to edit setup.py/xml files for your GPU/CUDA SDK
- Unify source: there is only one folder with all the CUDA code for both python and MATLAB. 
- Detector rotation now supported on parallel beam
## BugFix

- Fix bug where compute 5.0 was not being compiled in python
- Some arbitrary rotation angles were not being passed to the kernel 
- Reduced memory needed in MATLAB NikonDataLoader