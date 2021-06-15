Changelog
======

Please, add here all changes added since previous version of TIGRE (v2.1).
Many contributors have made this possible. Sincere thanks.

## Features

- You can now select which GPUs to run the code on, either by name, or GPUid.
- Improved compilation, now there is no need to edit setup.py/xml files for your GPU/CUDA SDK
- Unify source: there is only one folder with all the CUDA code for both python and MATLAB. 
- Detector rotation now supported on parallel beam
- Add NikonDataLoader() to python. By outside github contrib.
- Code quality optional code added, to be able to maintain decent python.
- Demos for python available
- Add 3D Shepp-Logan head phantom to Python.

## BugFix

- Fix bug where compute 5.0 was not being compiled in python
- Some arbitrary rotation angles were not being passed to the kernel 
- Reduced memory needed in MATLAB NikonDataLoader
- Mayor bug that impedes many iterative algorithms to work in python. The way order_subsets returns arrays changed, now they are not dtype=object, but dtype=float
- ComputeV errored for triplet angle inputs, fixed
- Fixed many python algorithm bugs (lost track of some)
- Removed many unnecesary files
- Fix and refactor checkDevices
- Fix backprojection type variable name and its default value in Python.