# pCT extension

## Introduction  

Cuda code and Matlab top layer to create optimized proton radiographs based on the method of [Collins-Fekete](https://doi.org/10.1088/0031-9155/61/23/8232). 
The radiograpies obtained from the code extension are intended to be used as pre-processing step before using the reconstruction algorithms implemented in TIGRE to do a full pCT reconstruction.
However, the the code could be also used as stand-alone tool for optimized proton radiograpies only.

## Manual

### Start up procedure

To use the pCT expansion within the TIGRE framework the CUDA source files have to be compiled. 
Starting from the MATLAB folder within the TIGRE framework, run the following script:

``` 
Compile_pCT (parallel geometry)
Compile_pCT_cone (cone beam geometry)
```

You will be asked to enter a pixel size and/or a intercepts vector size. If you want to stick to the default values, simply hit enter. The code will give you feedback on how many protons were rejected during the calculation of the radiographies due to the intercepts vector size.

### Data I/O
The measurement/simualtion data corresponding to the pCT experiment can be be forwarded as any Matlab-supported data type (.mat/.raw are recommended).

The simulation data has to contain the following parameters

| Parameter | short name | unit |
|---        |---         |---   |
| water equivalent path length      | WEPL    | cm|
| entry position (x,y) component   | posIn   | mm/cm/..|
| exit position (x,y) component | posOut | mm/cm/..|
| entry direction (x,y) component   |   dirIn| mm/cm/.., dz = 1!| 
| exit direction (x,y) component   | dirOut| mm/cm/.., dz = 1!|

The WEPL has to be given in cm. All other inputs have to share the same unit of length, which can be freely chosen. However, dz has to be 1 in all cases.

### Main method 
To ensure that the method delivers reasonable results particles with large angular and energetic deviations have to be cut from the dataset (3 sigma cuts) beforehand.
The main method for generating the improved radiographies is called by, 

```matlab
proj = pCTCubicSpline_mex(posIn, posOut, dirIn, dirOut, Wepl, PIXEL_SIZE, DETECTOR_SIZE_X, DETECTOR_SIZE_Y, POS_DET_IN, POS_DET_OUT, E_INIT, HULL_PARAMS);
```

While __Wepl__ is a __1D array of size N__ containing the WEPL measurement for each proton, __posIn, posOut, dirIn, dirOut__ are __1D arrays of size 2*N__ where first the x entries of all protons are filled in and secondly all y entries.
Other parameters that have to be defined are

| Parameter | Description | Unit | Type |
|---        |---         |---   |---  |
| PIXEL_SIZE      | side length of pixels used in the optimized projection | mm/cm/.. | single |
| DETECTOR_SIZE_X   | x length of the projection (in mm/cm/..)   | mm/cm/.. | int |
| DETECTOR_SIZE_Y | y length of the projection (in mm/cm/..) | mm/cm/.. | int |
| POS_DET_IN | position of the entry detector (always < 0!) | mm/cm/..| single |
| POS_DET_OUT | position of the exit detector (always > 0!)| mm/cm/..| single |
| E_INIT | initital proton energy | MeV | single |
| HULL_PARAMS | array of convex hull parameters [a, b, alpha, h] | mm/cm/.., alpha in rad | single |

If parameter __h__ of the HULL_PARAMS is set to zero, to calculation will be performed without a convex hull. Otherwise, the hull follows the equation (x cos(alpha)-z sin(alpha))²/a² + (x sin(alpha)+z cos(alpha))²/b² = 1 for y in [-h/2, h/2]. The ellispe can be rotated around the y-axis by the angle alpha.
Please click ![here](pct_cuda_toolbox.pdf "pCT toolbox geometry") to find a sketch displaying the definition of the geometry parameters used within the code extension. The origin is located in the center of the phantom.
Furthermore, a sketch displaying the definition of the detector parameters used within the code extension can be found ![here](pct_cuda_toolbox2.pdf "pCT toolbox detector parameters").

### Testing data
A demo script was prepared for the parallel geometry. To check if the code runs properly on your system run the script pCT_demo.m in MATLAB/Demos. A circle should be displayed if everything works correctly.

## Debugging notes

- If the general compilation of the pCT code with ```Compile_pCT``` fails with the compiler hinting at missing symbols, check if libstdc++ shipped with MATLAB is compatible with the CUDA drivers and the systems gcc.


Code tested and executed with Debian 10, CUDA 11.2, GPU Nvidia RTX 4000 and MATLAB R2020b and R2021b.

A python header is work-in-progress.


