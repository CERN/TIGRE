# Python

TIGRE supports Python 3.11 to 3.14 on Windows and Linux. We recommend
installing it in its own virtual environment (venv or conda).

## Prebuilt wheels (Linux)

Each [GitHub release](https://github.com/CERN/TIGRE/releases) has prebuilt
wheels for x86-64 Linux. They are built with CUDA 12.8 and bundle the CUDA
runtime, so you need only an NVIDIA driver that supports CUDA 12.8; no CUDA
toolkit or compiler.

1. Open the [latest release](https://github.com/CERN/TIGRE/releases/latest)
   and download the wheel that matches your Python version (`cp312` for
   Python 3.12, and so on).
2. Install it:

   ```sh
   pip install tigre-X.Y.Z-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
   ```

If there is no wheel for your Python version, install from source instead.

## Install from source

Installing from source compiles TIGRE's CUDA code on your machine, which takes
a few minutes.

### Requirements

::::::{tab-set}

::::{tab-item} Windows
- Python 3.11–3.14
- An NVIDIA GPU, or an AMD GPU with ROCm (see [AMD GPUs](#amd-gpus-rocm))
- The [CUDA toolkit](https://developer.nvidia.com/cuda-downloads)
- [Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/)
  with the C++ workload and the Windows SDK

Tested on Windows 10 and 11.

:::{note}
Your Windows user name must not contain spaces.
:::
::::

::::{tab-item} Linux
- Python 3.11–3.14
- An NVIDIA GPU, or an AMD GPU with ROCm (see [AMD GPUs](#amd-gpus-rocm))
- The [CUDA toolkit](https://developer.nvidia.com/cuda-downloads); follow
  [NVIDIA's Linux installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/)
  carefully
- gcc, in a version your CUDA toolkit supports

Tested on Ubuntu 16.04 and newer.
::::
::::::

### Straight from GitHub

```sh
# latest master
pip install git+https://github.com/CERN/TIGRE.git

# a specific release
pip install git+https://github.com/CERN/TIGRE.git@vX.Y.Z
```

Replace `vX.Y.Z` with a release tag. The latest release is **v{{ release }}**.

### From a local clone

Use this if you want to change TIGRE's code, or run the demos and tests.

```sh
git clone https://github.com/CERN/TIGRE.git
cd TIGRE
pip install .
```

Use `pip install -e .` for an editable install while developing. See
[contributing.md](https://github.com/CERN/TIGRE/blob/master/contributing.md)
for the development setup.

### If TIGRE cannot find CUDA

The build looks for CUDA in this order: the `CUDA_HOME` environment variable,
then `CUDA_PATH`, then `nvcc` on your `PATH`, then the default install location
(`C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\vX.Y` on Windows,
`/usr/local/cuda` on Linux). If it fails, point it at your installation:

::::{tab-set}
:::{tab-item} Windows
```bat
set CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8
pip install .
```
:::
:::{tab-item} Linux
```sh
export CUDA_HOME=/usr/local/cuda
pip install .
```
:::
::::

### Choosing the GPU backend

TIGRE is built for one GPU backend. When you build from source, it uses the
first one whose compiler it finds, in this order:

| Backend | GPUs | Found from |
|---|---|---|
| `cuda` | NVIDIA | `CUDA_HOME` or `CUDA_PATH`, `nvcc` on the `PATH`, or the default CUDA install location |
| `hip` | AMD | `ROCM_PATH` or `HIP_PATH`, `hipcc` on the `PATH`, or `/opt/rocm` |
| `sycl` | Intel | `icpx` (`icx` on Windows) on the `PATH` |
| `metal` | Apple Silicon | a Mac with Apple Silicon |

The build prints which backend it chose and why. To choose one yourself, set
`TIGRE_BACKEND`, for example `TIGRE_BACKEND=cuda pip install .`
(`set TIGRE_BACKEND=cuda` first on Windows).

At the moment `cuda` and `hip` can be built. The others are being added:
`sycl` in [#885](https://github.com/CERN/TIGRE/issues/885) and `metal` in
[#876](https://github.com/CERN/TIGRE/issues/876). If the build picks one of
them, it stops and says so; set `TIGRE_BACKEND=cuda` to build for NVIDIA GPUs.

### AMD GPUs (ROCm)

The `hip` backend builds TIGRE for Python on AMD GPUs (the MATLAB build is
CUDA only). It compiles the same `.cu` sources with `hipcc` from
[ROCm](https://rocm.docs.amd.com/) and links the HIP runtime (`amdhip64`)
instead of `cudart`. You need ROCm 7.2 or newer and a supported AMD GPU.

```sh
export TIGRE_BACKEND=hip
export ROCM_PATH=/opt/rocm
export HIP_ARCH=gfx90a
pip install .
```

- `ROCM_PATH` points to the ROCm install (default `/opt/rocm`, or the prefix
  of the `hipcc` on your `PATH`).
- `HIP_ARCH` is the target GPU architecture for `--offload-arch` (default
  `gfx90a`; the build prints the one it uses). Pass a comma-separated list to
  build for several, for example `HIP_ARCH=gfx90a,gfx1100,gfx1201`.
- On Windows, also set `HIPCC` to the ROCm `clang++` (for example
  `%ROCM_PATH%\lib\llvm\bin\clang++.exe`).

Tested on Linux with gfx90a (Instinct MI250X) and gfx1100 (Radeon PRO W7800),
and on Windows with gfx1201 (RDNA4) and gfx1151 (Radeon 8060S).

## Conda (unofficial)

CCPi Tomographic Imaging builds and maintains a conda package:

```sh
conda install -c ccpi tigre
```

It is not maintained by the TIGRE team, so please report problems with it to
the [TIGRE-conda repository](https://github.com/TomographicImaging/TIGRE-conda/issues).

## Check your installation

This lists your GPUs and runs a small reconstruction:

```python
import numpy as np
import tigre
import tigre.algorithms as algs
from tigre.utilities import gpu, sample_loader

print(gpu.getGpuNames())

geo = tigre.geometry_default(high_resolution=False)
angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
head = sample_loader.load_head_phantom(geo.nVoxel)
projections = tigre.Ax(head, geo, angles)
reconstruction = algs.fdk(projections, geo, angles)
print(reconstruction.shape)
```

If that works, continue with the demos in `Python/demos`, starting with
`d01_CreateGeometry.py`.

## Advanced: disable pinned memory

TIGRE uses pinned (page-locked) host memory to overlap memory transfers with
computation. If you reconstruct datasets larger than your RAM and rely on swap,
compile without pinned memory. The operators will be slower.

Set the `TIGRE_NO_PINNED_MEMORY` environment variable when you install from
source:

::::{tab-set}
:::{tab-item} Windows
```bat
set TIGRE_NO_PINNED_MEMORY=1
pip install --no-cache-dir git+https://github.com/CERN/TIGRE.git
```
:::
:::{tab-item} Linux
```sh
TIGRE_NO_PINNED_MEMORY=1 pip install --no-cache-dir git+https://github.com/CERN/TIGRE.git
```
:::
::::

`--no-cache-dir` stops pip from reusing a wheel it built earlier with pinned
memory. The same works with `pip install .` from a clone. The prebuilt wheels
always use pinned memory.
