Installation Instructions for Python
======

## Windows

### Requirements:

1. Python 3
2. MSVC
3. A CUDA capable GPU from NVIDIA with [compute capability greater or equal to 3.0](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
4. CUDA Toolkit

Tested on

| Software        | Version           | 
| ------------- |:-------------:|
|**Windows**| 10 |
|**Python**| 3.7 3.8 |
|**CUDA**| 10.1 |
|**MSVC**| 19.24 |

### Simple Instructions

1. Install Python and pip, MSVC and CUDA
2. run `git clone https://github.com/CERN/TIGRE.git` 
3. run `python setup.py install --user` in the pyhton folder. 

A succesfull installation should be able to execute the script at TIGRE/Python/example.py

###  Step by Step Instructions:

1. Install [MS Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/) with Windows SDK.

2. Install Python 3 + pip. Use a virtual [conda environment](https://www.anaconda.com/) or jus tnormal python installation.

3. Install [CUDA](https://developer.nvidia.com/cuda-downloads). Make sure the `CUDA_PATH` and
      `PATH` environment variable are set accordingly.
   
4. Download/clone TIGRE
    
	`git clone https://github.com/CERN/TIGRE.git` 

5. Compile libraries

	`cd TIGRE/Python/` 
	`python setup.py install --user`

	Install in this case will make a copy of pytigre to your python distribution. Therefore the `develop` command is more useful when modifying the source files and developing the software. 

	`python setup.py develop --user`

6. Try demo 3. If it runs succesfully then you are good to go.

**Note:** It is known that the package cannot be imported using a pure Python 3.8 installation, i.e.
not using a conda environment. Please try a pure installation of Python 3.7 or consider using a
conda environment.

## Linux

### Requirements:

1. Python 2/Python 3
2. gcc
3. A CUDA capable GPU from NVIDIA with [compute capability greater or equal to 3.0](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
4. CUDA Toolkit


Tested on

| Software        | Version           | 
| ------------- |:-------------:|
|**Ubuntu**| 16.04 17.10|
|**Python**| 2.7 3.7 |
|**CUDA**| 8.0 9.2 10.1 10.2|
|**gcc**|  7.6.0|

### Simple Instructions

1. Install python, gcc, pip and CUDA
2. run `git clone https://github.com/CERN/TIGRE.git` 
3. run `python setup.py install --user` in the pyhton folder. 

A succesfull installation should be able to execute the script at TIGRE/Python/example.py


###  Step by Step Instructions:

For Ubuntu

1. Install python  and pip (you can use 2 or 3, example show for 2)

	```
	sudo apt update
	sudo apt upgrade
	sudo apt install python2.7 python-pip
	```
	
2. Install CUDA

   Installing CUDA in linux (specially one with a GUI) can be a challenge. Please follow [NVIDIAs instructions](https://developer.download.nvidia.com/compute/cuda/10.0/Prod/docs/sidebar/CUDA_Installation_Guide_Linux.pdf) carefully.\
   [CUDA download link](https://developer.nvidia.com/cuda-downloads)

3. Install gcc 

   gcc shoudl already be installed in your linux, as it is part of the linux distribution.\
   If you need to install an older version of gcc, [read here](https://askubuntu.com/questions/923337/installing-an-older-gcc-version3-4-3-on-ubuntu-14-04-currently-4-8-installed).
   
4. Download/clone TIGRE
    
	`git clone https://github.com/CERN/TIGRE.git` 

5. Compile libraries

	`cd TIGRE/Python/` 
	`python setup.py install --user`

	Install in this case will make a copy of pytigre to your python distribution. Therefore the `develop` command is more useful when modifying the source files and developing the software. 

	`python setup.py develop --user`

6. Try demo 3. If it runs succesfully then you are good to go. 



if this fails, then try:
	
`export CUDAHOME=yourcudahome`, e.g. default is `export CUDAHOME=/usr/local/cuda`
`python setup.py install --user`
	
	
**NOTE:** as of November 2020 the pip pytigre is behind the main repo, we recomedn you install it and compile it yourself. Trying to fix that. 