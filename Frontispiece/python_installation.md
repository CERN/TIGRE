# Installation Instructions for Python

## Table of contents

- [Windows](#windows)
   - [Requirements](#requirements)
   - [Build Instructions](#build-instructions)
   - [Conda package](#conda-install)
- [Linux](#linux)
   - [Requirements](#requirements-1)
   - [Build Instructions](#build-instructions-1)
   - [Conda package](#conda-install-1)
- [Optional Code Style Enforcement](#optional-code-style-enforcement)
- [Advanced](#advanced)

## Windows

### Requirements

1. Python 3.7-3.11
2. NVIDIA CUDA-capable GPU with [compute capability >=3.5](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
3. MSVC (for building)
4. CUDA Toolkit (for building)

Tested on:

| Software      | Version       | 
| ------------- |:-------------:|
|**Windows**| 10, 11 |
|**Python**| 3.10-3.14 |
|**CUDA**| >=9.2 |
|**MSVC**| 19.24 |

### Build Instructions

We strongly recommend using `conda` environments and doing the install in one specific for tigre. 

> [!TIP]
> Quick Summary:
>
> 1. Install Python and pip, MSVC and CUDA
> 2. run `git clone https://github.com/CERN/TIGRE.git` 
> 3. run `pip install .` in the root folder. 

A successful installation should be able to execute the script at `TIGRE/Python/example.py`

1. Install [MS Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/) with Windows SDK.\
   **NOTE:** The User has to have no spaces.

2. Install Python 3 + pip. You can use a virtual [conda environment](https://www.anaconda.com/) or just a normal python installation. We strongly recommend the conda environments. Activate the environment before step 5.

3. Install [CUDA](https://developer.nvidia.com/cuda-downloads). Make sure the `CUDA_PATH` and `PATH` environment variable are set accordingly.\
  **NOTE:** The User has to have no spaces.

4. Download/clone TIGRE

	`git clone https://github.com/CERN/TIGRE.git` 

5. Compile libraries

	```sh
	cd TIGRE/  
	pip install . --user
	```
	**NOTE:** If you are working under the virtual environment that created by `venv` or a `conda` environment and you want to install TIGRE to it, 
	you should remove the `--user` option. 
	With the `--user` option, TIGRE and the other required packages will be installed to your Python user install directory, not to your virtual environment or system directory.

6. Try demo 3. If it runs successfully then you are good to go.

### Conda install

An unofficial conda package is built & maintained by CCPi Tomographic Imaging - you can install it by running:

`conda install -c ccpi tigre`

Please report conda issues [upstream](https://github.com/TomographicImaging/TIGRE-conda/issues).

## Linux

### Requirements

We strongly recommend using `conda` environments and doing the install in one specific for tigre. 

1. Python 3
2. NVIDIA CUDA-capable GPU with [compute capability >=3.5](https://en.wikipedia.org/wiki/CUDA#GPUs_supported)
3. gcc (for building)
4. CUDA Toolkit (for building)

Tested on:

| Software      | Version       | 
| ------------- |:-------------:|
|**Ubuntu**| >=16.04 |
|**Python**| 3.10-3.14 |
|**CUDA**| >=9.2 |
|**gcc**| 7.6.0 |

### Build Instructions

We strongly recommend using `conda` environments and doing the install in one specific for tigre. 

> [!TIP]
> Quick Summary:
>
> 1. Install python and pip, gcc and CUDA
> 2. run `git clone https://github.com/CERN/TIGRE.git` 
> 3. run `pip install .` in the root folder. 

A successful installation should be able to execute the script at `TIGRE/Python/example.py`

For Ubuntu

1. Install python and pip

	Recommended to do it via [Anaconda3](https://docs.anaconda.com/free/anaconda/install/linux/), and not the following. 

	```sh
	sudo apt update
	sudo apt upgrade
	sudo apt install python3.10 python-pip
	```

2. Install CUDA

   Installing CUDA in linux (specially one with a GUI) can be a challenge. Please follow [NVIDIAs instructions](https://developer.download.nvidia.com/compute/cuda/10.0/Prod/docs/sidebar/CUDA_Installation_Guide_Linux.pdf) carefully.\
   [CUDA download link](https://developer.nvidia.com/cuda-downloads)

3. Install gcc

   gcc should already be installed in your linux, as it is part of the linux distribution.\
   If you need to install an older version of gcc, [read here](https://askubuntu.com/questions/923337/installing-an-older-gcc-version3-4-3-on-ubuntu-14-04-currently-4-8-installed).

4. Download/clone TIGRE

	`git clone https://github.com/CERN/TIGRE.git` 

5. Compile libraries

	```
	cd TIGRE/  
	pip install . --user
	```
	**NOTE:** If you are working under the virtual environment that created by `venv` and you want to install TIGRE to it, 
	you should remove the `--user` option. 
	With the `--user` option, TIGRE and the other required packages will be installed to your Python user install directory, not to your virtual environment or system directory.

6. Try demo 3. If it runs successfully then you are good to go. 

if this fails, then try:

`export CUDAHOME=yourcudahome`, e.g. default is `export CUDAHOME=/usr/local/cuda`
`pip install . --user`

**NOTE:** as of November 2020 the pip pytigre is behind the main repo, we recommend you install it and compile it yourself. Trying to fix that. 

### Conda install

An unofficial conda package is built & maintained by CCPi Tomographic Imaging - you can install it by running:

`conda install -c ccpi tigre`

Please report conda issues [upstream](https://github.com/TomographicImaging/TIGRE-conda/issues).

## Optional Code Style Enforcement
Optional linting dependencies are provided to enforce the prevailing codestyle in the Python component of the TIGRE library.

The primary linting packages utilized are:

* [black](https://black.readthedocs.io/en/stable/) for automatic code formatting
* [flake8](https://flake8.pycqa.org/en/latest/) for code style enforcement

Multiple installation options have been provided:

### `pip` Installation
A `dev` dependency group is included in the `pyproject.toml` file. Use `pip install --group dev .` to install those dependencies. 

### pre-commit
Linting dependencies are also provided as a [pre-commit](https://pre-commit.com/) configuration. With the pre-commit hooks installed, code will be linted prior to commit, preventing code from being committed that fails linting. By default, these hooks are only run against staged files.

The pre-commit package is included in the above `pip` installs; once the framework is installed, hooks can be installed by invoking `pre-commit install`.

### Usage
Once installed, invoking the linting tools is straightforward. From TIGRE's `Python` directory:

```bash
# Run in current directory & subdirectories
$ flake8 .
$ black .

# Or specify a file
$ flake8 ./example.py
$ black ./example.py
```

pre-commit may also be manually invoked, which will run `flake8` and `black` against all staged files:

```bash
$ pre-commit run
flake8...................................................................Passed
black....................................................................Passed
Check for merge conflicts................................................Passed
Check Toml...........................................(no files to check)Skipped
Check Yaml...........................................(no files to check)Skipped
Fix End of Files.........................................................Passed
Mixed line ending........................................................Passed
check blanket noqa.......................................................Passed
```

**NOTE:** pre-commit may also be manually invoked against *all* files (staged and unstaged) using the `pre-commit run --all-files`. However, some changes made to Python's TIGRE codebase by `black` have been manually reverted for readability reasons and should not be committed in their blackened state.

## Advanced 

If you are doing reconstruction of large datasets, and you want to use swap memory, you will need to deactivate TIGREs pinned memory feature at compile time. This will allow you to use swap memory, but it will make the operators in TIGRE slower, as pinned memory is used for simultaneous memory and compute. 

You can do this by calling the `setup.py` with the flag `--no_pinned_memory`. 
