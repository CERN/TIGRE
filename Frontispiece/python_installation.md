Installation Instructions for Python
======

## Windows

The python version of TIGRE has not been teste on a windows machine yet. Please feel free to test it and help us build an easy installer for this OS. [Contact us](mailto:ander.biguri@gmail.com).

## Linux

### Requirements:

1. Python 2.7
2. gcc
3. A CUDA capable GPU from NVIDIA with [compute capability greater or equal to 3.0](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
4. CUDA Toolkit


Tested on

| Software        | Version           | 
| ------------- |:-------------:|
|**Ubuntu**| 16.04 17.10|
|**Python**| 2.7 |
|**CUDA**| 8.0 9.2 |
|**gcc**|  7.6.0|

### Simple Instructions

1. Install python, gcc, pip and CUDA
2. run `pip install pytigre`
3. Run `python setup.py install`. 

A succesfull installation should be able to execute the script at  XXX


###  Step by Step Instructions:

For Ubuntu

1. Install python 2.7 and pip

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
   
4. Get TIGRE
    
	`pip install pytigre` 
	
6. Try demo 3. If it runs succesfully then you are good to go.

Instead, if you rather compile from source, download/clone the repository and then run 

 
	`python setup.py install` 
	
	if this fails, then try:
	
	`export CUDAHOME=yourcudahome`, e.g. default is `export CUDAHOME=/usr/local/cuda`
	`sudo --preserve-env python setup.py install`
	
	
