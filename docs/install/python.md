---
og:description: "Install TIGRE for Python: prebuilt wheels for Linux with one pip command, or building from source on Windows and Linux."
myst:
  html_meta:
    description: "Install TIGRE for Python: prebuilt wheels for Linux with one pip command, or building from source on Windows and Linux."
---

# Python

TIGRE supports Python 3.11 to 3.14 on Windows and Linux. We recommend
installing it in its own virtual environment (venv or conda).

## Prebuilt wheels (Linux)

Each [GitHub release](https://github.com/CERN/TIGRE/releases) has prebuilt
wheels for x86-64 Linux. They bundle the CUDA runtime, so you need no CUDA
toolkit or compiler: only an NVIDIA GPU and its driver (see
[what the wheels need](#what-the-wheels-need)).

To install TIGRE {{ latest_version }}, the latest release, pick your Python
version and run the command in your environment (`python --version` tells you
which you have):

```{include} ../_generated/wheels.md
```

Then [check that it works](#check-your-installation). If there is no wheel for
your Python version, or you are on Windows, install from source instead, or try
the [conda package](#conda-unofficial).

### What the wheels need

- **Linux on x86-64** with glibc 2.28 or newer: Ubuntu 18.10+, Debian 10+,
  RHEL/Rocky/Alma 8+, Fedora 29+ and most current distributions.
- **An NVIDIA GPU of compute capability 5.0 or newer** (Maxwell, from 2014,
  to Blackwell). Look your GPU up on
  [NVIDIA's list](https://developer.nvidia.com/cuda-gpus).
- **An NVIDIA driver for CUDA 12.8**: version 570 or newer. Drivers from
  525 on support all of CUDA 12 through NVIDIA's
  [minor-version compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/)
  and should also work, but have not been tested with TIGRE.
  `nvidia-smi` prints your driver version.

## Install from source

Installing from source compiles TIGRE's CUDA code on your machine, which takes
a few minutes.

### Requirements

::::::{tab-set}

::::{tab-item} Windows
- Python 3.11–3.14
- An NVIDIA GPU
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
- An NVIDIA GPU
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

At the moment only `cuda` can be built. The others are being added: `hip` in
[#747](https://github.com/CERN/TIGRE/pull/747), `sycl` in
[#885](https://github.com/CERN/TIGRE/issues/885) and `metal` in
[#876](https://github.com/CERN/TIGRE/issues/876). If the build picks one of
them, it stops and says so; set `TIGRE_BACKEND=cuda` to build for NVIDIA GPUs.

## Conda (unofficial)

CCPi Tomographic Imaging builds and maintains a conda package:

```sh
conda install -c ccpi tigre
```

It is not maintained by the TIGRE team, so please report problems with it to
the [TIGRE-conda repository](https://github.com/TomographicImaging/TIGRE-conda/issues).

## Check your installation

This lists your GPUs and runs a small reconstruction. Run it from a folder
outside TIGRE's source code: inside it, Python finds the uncompiled sources
instead of the installed package.

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

It should print your GPUs' names and then `(64, 64, 64)`. If
that works, continue with the [demos](../demos/index.rst). If it fails, see
[](../user-guide/troubleshooting.md#installing).

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
