Installation Instructions for MATLAB
======

## Table of contents

- [Windows](#windows)
   - [Requirements](#requirements)
   - [Simple Instructions](#simple-instructions)
   - [Step by Step Instructions](#step-by-step-instructions)
- [Linux](#linux)
   - [Requirements](#requirements-1)
   - [Simple Instructions](#simple-instructions-1)
   - [Step by Step Instructions](#step-by-step-instructions1)
- [Advanced](#advanced)
   *****
   
## Windows

### Requirements:

1. MATLAB
2. Visual Studio (Community or Profesional)
3. A CUDA capable GPU from NVIDIA with [compute capability greater or equal to 3.5](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
4. CUDA Toolkit (9.2 or newer)

Tested on

| Software        | Version           | 
| ------------- |:-------------:|
|**Windows**| 7, 8, 10.|
|**MATLAB**|Any MATLAB >2016b|
|**CUDA**|Any CUDA 9.2>|
|**Visual Studio**| 2010 2013 2015 2019 2022|



### Simple Instructions

1. Install MATLAB, Visual Studio and CUDA (Remember to install C++ when isntalling Visual Studio!)
2. Run `Compile.m`

A succesfull installation should be able to execute the script at `TIGRE/MATLAB/Demos/d03_generateData.m` without errors.



###  Step by Step Instructions:

1. Install MATLAB

2. Install CUDA

   Any version avobe 9.2 has been tested, however its recommended to get the latests version as possible, for performance and support.\
   [CUDA download link](https://developer.nvidia.com/cuda-downloads)\
   [Detailed installation guide](https://developer.download.nvidia.com/compute/cuda/10.0/Prod/docs/sidebar/CUDA_Installation_Guide_Windows.pdf)\
   **NOTE**: In windows at least, the User has to have no spaces. 
   
3. Install Visiual Studio

   **Make sure you install C++**.\
   [Download link for the latest version](https://visualstudio.microsoft.com/downloads/)\
   [Download link for older versions](https://visualstudio.microsoft.com/vs/older-downloads/)\
   **NOTE**: In windows at least, the User has to have no spaces. 
   
4. Download TIGRE

   If you are using git, run: `git clone https://github.com/CERN/TIGRE.git`\
   Manually [download zip file](https://github.com/CERN/TIGRE/archive/master.zip) oherwise.

5. Test the correct configuration of Visual Studio 

   Open MATLAB and run `mex -setup -v`. Among other things, the output should contain:
   ```
   ... Looking for compiler 'Microsoft Visual C++ 2015 (C)' ...
   ... Looking for registry setting 'HKLM\SOFTWARE\Microsoft\VisualStudio\SxS\VC7' 14.0 ...No.
   ... Looking for registry setting 'HKCU\SOFTWARE\Microsoft\VisualStudio\SxS\VC7' 14.0 ...No.
   ... Looking for registry setting 'HKLM\SOFTWARE\Wow6432Node\Microsoft\VisualStudio\SxS\VC7' 14.0 ...Yes ('C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\').
   ... Looking for file 'C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\cl.exe' ...Yes.
   ... Looking for folder 'C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC' ...Yes.
   ```

   Note the last 3 **Yes**. If those are not there it means that you do not have installed C++ on step 3 of this tutorial.\
   To fix:\
   -go to `Control panel>Add or remove programs -> Visual studio community 20XX -> modify`\
   -`Lanugages -> Visual C++`
   
   Make sure that when you run `mex -setup -v` C++ is installed and Visual Studio is selected as the compiler for C/C++
   
6. By opening MATLAB on `yourTIGREpath/MATLAB`, execute `Compile.m`

   If it fails, try opening `mex_CUDA_win64.xml` with your favourite editor and changing line 125 to link to your local `nvcc` location.
   
7. Initialize TIGRE by typing `InitTIGRE` on the MATLAB Command Window.

8. Run file `TIGRE/MATLAB/Demos/d03_generateData.m`. If it successfully executes, you have installed and compiled TIGRE properly.


If none of this works, please contact the authors at [tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com) or [ander.biguri@gmail.com](mailto:ander.biguri@gmail.com)
   
   
## Linux

### Requirements:
1. MATLAB
2. gcc
3. A CUDA capable GPU from NVIDIA with [compute capability greater or equal to 3.0](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
4. CUDA Toolkit (9.2 or newer)

Tested on


| Software       | Version    | 
| ------------- |:-------------:|
| **Ubuntu**|  Any ubuntu 16.04>| 
| **MATLAB**|  Any MATLAB 2016b>| 
| **CUDA**| Any Cuda 0.2>| 
| **gcc**|  6.4.0 7.2.0| 

### Simple Instructions

1. Install MATLAB, gcc and CUDA
2. Run `Compile.m`

A succesfull installation should be able to execute the script at `TIGRE/MATLAB/Demos/d03_generateData.m` without errors.

### Step by Step Instructions:

1. Install MATLAB

2. Install CUDA

   Installing CUDA in linux (specially one with a GUI) can be a challenge. Please follow [NVIDIAs instructions](https://developer.download.nvidia.com/compute/cuda/10.0/Prod/docs/sidebar/CUDA_Installation_Guide_Linux.pdf) carefully.\
   [CUDA download link](https://developer.nvidia.com/cuda-downloads)

3. Install gcc 

   gcc should already be installed in your linux, as it is part of the linux distribution.\
   If you need to install an older version of gcc, [read here](https://askubuntu.com/questions/923337/installing-an-older-gcc-version3-4-3-on-ubuntu-14-04-currently-4-8-installed).
   
4. Download TIGRE

   If you are using git, run: `git clone https://github.com/CERN/TIGRE.git`\
   Manually [download zip file](https://github.com/CERN/TIGRE/archive/master.zip) oherwise.
   
5. Make sure your terminal knows where CUDA is.

   For Ubuntu:
   - Using your favourite test editor, open ~/.bashrc. e.g. `gedit ~/.bashrc`.
   - Append to the file the following lines:
     ```
	 export LD_LIBRARY_PATH=/usr/local/cuda/lib64
     export PATH=$PATH:/usr/local/cuda/bin
	 ```
     
   - restart your terminal
   
6. From a terminal, execute MATLAB as `matlab` and run `Compile.m` located on `TIGRE/MATLAB` 

7. Run file `TIGRE/MATLAB/Demos/d03_generateData.m`. If it successfully executes, you have installed and compiled TIGRE properly.


If none of this works, please contact the authors at [tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com) or [ander.biguri@gmail.com](mailto:ander.biguri@gmail.com)

****

## Advanced

If you are doing reconstruction of large datasets, and you want to use swap memory, you will need to deactivate TIGREs pinned memory feature at compile time. This will allow you to use swap memory, but it will make the operators in TIGRE slower, as pinned memory is used for simultaneous memory and compute. 


You can do this by calling the `Compile.m` file from the MATLAB command line as `Compile --no_pinned_memory`.




