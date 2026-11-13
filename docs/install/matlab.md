# MATLAB

TIGRE for MATLAB is compiled on your machine with `Compile.m`.

## Requirements

::::::{tab-set}

::::{tab-item} Windows
- MATLAB R2016b or newer
- An NVIDIA GPU
- The [CUDA toolkit](https://developer.nvidia.com/cuda-downloads)
- Visual Studio 2013, 2015, 2017, 2019 or 2022, **with the C++ workload**
  ([latest](https://visualstudio.microsoft.com/downloads/),
  [older versions](https://visualstudio.microsoft.com/vs/older-downloads/))

:::{note}
Your Windows user name must not contain spaces.
:::
::::

::::{tab-item} Linux
- MATLAB R2016b or newer
- An NVIDIA GPU
- The [CUDA toolkit](https://developer.nvidia.com/cuda-downloads); follow
  [NVIDIA's Linux installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/)
  carefully
- gcc, in a version supported by both your MATLAB release and your CUDA toolkit
::::
::::::

## Install

1. Get TIGRE, with `git clone https://github.com/CERN/TIGRE.git` or by
   downloading the [zip file](https://github.com/CERN/TIGRE/archive/master.zip).

2. **Linux only:** make sure MATLAB can find CUDA. Add these lines to
   `~/.bashrc`, open a new terminal, and start MATLAB from it:

   ```sh
   export PATH=$PATH:/usr/local/cuda/bin
   export LD_LIBRARY_PATH=/usr/local/cuda/lib64
   ```

3. **Windows only:** check that MATLAB sees Visual Studio's C++ compiler.
   Run `mex -setup -v` in MATLAB. Near the end, the lines that look for
   Visual C++ should end in **Yes**, and Visual Studio should be selected as
   the C/C++ compiler. If not, open the Visual Studio Installer, choose
   **Modify**, and add the C++ workload.

4. In MATLAB, go to the `TIGRE/MATLAB` folder and run:

   ```matlab
   Compile
   ```

   On Windows, `Compile` picks the configuration file matching your Visual
   Studio version (`mex_CUDA_win64_MVS20XX.xml`). To force a specific one, run
   `Compile -mex_xml mex_CUDA_win64_MVS2019.xml`.

5. Add TIGRE to the MATLAB path. From the `TIGRE/MATLAB` folder, run:

   ```matlab
   InitTIGRE
   ```

   `InitTIGRE` adds folders relative to the current folder, so run it from
   `TIGRE/MATLAB` at the start of every MATLAB session that uses TIGRE.

6. Check the installation by running `Demos/d03_generateData.m`. If it runs
   without errors, TIGRE is installed.

## Advanced: disable pinned memory

TIGRE uses pinned (page-locked) host memory to overlap memory transfers with
computation. If you reconstruct datasets larger than your RAM and rely on swap,
compile without pinned memory. The operators will be slower.

```matlab
Compile -no_pinned_mem
```

## Getting help

If installation fails, search the
[GitHub issues](https://github.com/CERN/TIGRE/issues?q=is%3Aissue) and
[Discussions](https://github.com/CERN/TIGRE/discussions), then ask there,
including your OS, GPU, CUDA version, compiler version and MATLAB release.
You can also email [tigre.toolbox@gmail.com](mailto:tigre.toolbox@gmail.com).
