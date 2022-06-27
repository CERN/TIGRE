# pCT extension

## Introduction  

Cuda code and Matlab top layer to create optimized proton radiographies based on the method of [Collins-Fekete](https://doi.org/10.1088/0031-9155/61/23/8232). 
The radiograpies obtained from the code extension are intended to be used as pre-processing step before using the reconstruction algorithms implemented in TIGRE to do a full pCT reconstruction.
However, the the code could be also used as stand-alone tool for optimized proton radiograpies only.

## Manual

### Start up procedure

To use the pCT extension within the TIGRE framework, the CUDA source files have to be compiled. 
Starting from the MATLAB folder within the TIGRE framework, run the following script:

``` 
CompilePCT 
```

You will be asked to enter a pixel length and/or a intercepts vector size. In case X/Y pixel lengths are not the same, enter the smaller one. If you want to stick to the default values, simply hit enter. The code will give you feedback on how many protons were rejected during the calculation of the radiographies due to the intercepts vector size.

### Data I/O
The measurement/simualtion data corresponding to the pCT experiment can be be forwarded as any Matlab-supported data type (.mat/.raw are recommended).

The simulation data has to contain the following parameters

| Parameter | short name | unit |
|---        |---         |---   |
| water equivalent path length      | WEPL    | mm |
| entry position (x,y) component   | posIn   | mm |
| exit position (x,y) component | posOut | mm |
| entry direction (x,y) component   |   dirIn| mm , dz = 1!| 
| exit direction (x,y) component   | dirOut| mm, dz = 1!|

All inputs have to share the same unit of length, which is mm. Be aware that dz has to be 1 in all cases.

### Main method 
To ensure that the method delivers reasonable results, particles with large angular and energetic deviations have to be cut from the dataset (3 sigma cuts) beforehand.
The main method for generating the improved radiographies is called by, 

```matlab
proj = pCTCubicSpline_mex(posIn, posOut, dirIn, dirOut, Wepl, eIn, geo);
```

While __Wepl__ is a __1D array of size N__ containing the WEPL measurement for each proton, __posIn, posOut, dirIn, dirOut__ are __1D arrays of size 2*N__ where first the x entries of all protons are filled in and secondly all y entries.
__eIn__ is a single value containing the initial beam energy (in MeV). Other parameters that have to be defined are following the geometry definition in TIGRE with additional parameters

| Parameter | Description | Unit | Type |
|---        |---         |---   |---  |
| geo.dDetector      | side length of pixels used in the optimized projection | mm | single |
| geo.sDetecor   | x and y length of the projection (in mm/cm/..) (array with two entries)  | mm | single |
| geo.DSO | distance between source and origin (always > 0!) | mm | single |
| geo.DSID | distance between source and entry detector (always > 0!) | mm | single |
| geo.DSD | distance between source and exit detector (always > 0!)| mm | single |
| geo.mode | "parallel" or "cone"  | / | string |
| geo.hull | array of convex hull parameters [a, b, alpha, h] | mm | single |

If parameter __h__ of the HULL_PARAMS is set to zero, to calculation will be performed without a convex hull. Otherwise, the hull follows the equation (x cos(alpha)-z sin(alpha))²/a² + (x sin(alpha)+z cos(alpha))²/b² = 1 for y in [-h/2, h/2]. The ellispe can be rotated around the y-axis by the angle alpha.
Please click ![here](pct_cuda_toolbox.pdf "pCT toolbox geometry") to find a sketch displaying the definition of the geometry parameters used within the code extension. The origin is located in the center of the phantom.
Furthermore, a sketch displaying the definition of the detector parameters used within the code extension can be found ![here](pct_cuda_toolbox2.pdf "pCT toolbox detector parameters").

### Testing data
A demo script calculating one optimized radiography was prepared for the parallel geometry. To check if the code runs properly on your system run the script pRad_demo.m in MATLAB/Demos. 

## Debugging notes

- If the general compilation of the pCT code with ```CompilePCT``` fails with the compiler hinting at missing symbols, check if libstdc++ shipped with MATLAB is compatible with the CUDA drivers and the systems gcc.


Code tested and executed with Debian 10, CUDA 11.2, GPU Nvidia RTX 4000 and MATLAB R2020b and R2021b.

A python header is work-in-progress.


