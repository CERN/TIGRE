---
og:description: "Fixes for common TIGRE problems: installation and compilation errors, GPU and memory errors, and frequently asked questions."
myst:
  html_meta:
    description: "Fixes for common TIGRE problems: installation and compilation errors, GPU and memory errors, and frequently asked questions."
---

# Troubleshooting and FAQ

If your problem is not covered here, search the
[issues](https://github.com/CERN/TIGRE/issues) and
[discussions](https://github.com/CERN/TIGRE/discussions). When you report a
bug, include a
[minimal, complete example](https://stackoverflow.com/help/minimal-reproducible-example)
that someone else can run to see the error, along with your operating system,
GPU, CUDA version and TIGRE version.

## Installing

**TIGRE does not compile.**
Follow the [installation guide](../install/index.md) for your language; it
covers the usual causes (CUDA not found, wrong compiler, missing Visual
Studio on Windows).

**`nvcc fatal: Unsupported gpu architecture 'compute_XX'`.**
Your CUDA version cannot build for one of the GPU architectures in the build
flags. The Python build asks `nvcc` which architectures it supports, and the
MATLAB build chooses them from your CUDA version, so this should be rare. If
you still see it, remove the unsupported `-gencode` entries from the
`COMPFLAGS` line of your MATLAB `mex_CUDA_*.xml` file, and please open an
issue so the build can be fixed.

**`ImportError: DLL load failed while importing _Ax` (Windows).**
TIGRE finds CUDA's DLLs through the `CUDA_PATH` environment variable. Set it
to the CUDA toolkit TIGRE was built with (for example
`C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8`), open a new
terminal, and try again. If several CUDA versions are installed, `CUDA_PATH`
must point to the one used for the build; otherwise reinstall TIGRE with it.

**`undefined symbol: cudaGetDeviceProperties_v2` or another `undefined symbol`
error from `_Ax` (Linux).**
Python loaded an older CUDA runtime than the one TIGRE was built with, often
from `LD_LIBRARY_PATH` or from a CUDA bundled in a conda environment. Point
`LD_LIBRARY_PATH` at the `lib64` folder of the CUDA used for the build, or
reinstall TIGRE in this environment so they match.

**An `ImportError` about `_Ax` or another compiled module, right after
installing.**
Check that you are not running Python from inside TIGRE's source folder (the
repository, or its `Python` folder): there, `import tigre` finds the
uncompiled sources instead of the installed package. Change to another folder
and try again.

**`CUDA driver version is insufficient for CUDA runtime version`.**
Your NVIDIA driver is older than the CUDA that TIGRE uses. Update the driver;
the prebuilt wheels need version 570 or newer (`nvidia-smi` shows yours).

**There are many warnings during compilation.**
They are expected and do not affect results or speed.

**Can I compile with MinGW on Windows?**
No. CUDA on Windows requires the Microsoft Visual Studio compiler. The free
Community edition is enough.

## Running

**`Input data should be float32` (Python) or `Image should be single type` (MATLAB).**
TIGRE does not convert arrays for you. Use `img.astype(np.float32)` or
`single(img)`.

**A shape error from `Ax`, `Atb` or an algorithm.**
The image must have shape `geo.nVoxel`, and the projections
`(len(angles), *geo.nDetector)` in Python. Remember that Python orders the
axes `(z, y, x)` and `(v, u)`; see [](geometry.md).

**`nVoxel*dVoxel is not equal to sVoxel`.**
When you change one of `nVoxel`, `dVoxel` or `sVoxel`, update the others to
match. The same applies to the detector fields.

**`texture memory fail`, `invalid device function` or `no kernel image is
available for execution on the device`.**
TIGRE was not built for your GPU. Rebuild it on the computer you use, with a
CUDA version that supports your GPU.

**`One (or more) of your GPUs is being heavily used by another program`.**
A warning: less than half of a selected GPU's memory is free, so TIGRE may
run out of it. Close programs using the GPU, or select another GPU; see
[](gpus.md). With [GPU arrays](interoperability.md#gpu-arrays) it is not
shown, as sharing the GPU (with PyTorch, for example) is expected there.

**On Windows, TIGRE or PyTorch becomes very slow when the GPU is nearly full.**
Recent NVIDIA drivers for Windows (536.40 and newer) move GPU memory to system
memory when the GPU runs out, instead of reporting "out of memory". Programs
keep running, but many times slower. To get the error instead, open the NVIDIA
Control Panel, go to *Manage 3D settings*, and set *CUDA - Sysmem Fallback
Policy* to *Prefer No Sysmem Fallback* (for all programs, or for your
`python.exe`). Then reduce what is on the GPU: the batch size, the image size,
or the number of angles per call.

**`the launch timed out and was terminated` or `unspecified launch failure`
on large images.**
The operating system stopped a GPU kernel because it ran too long, which
happens on GPUs that also drive a display. On Windows, raise the TDR delay or
turn TDR off (for example in NVIDIA Nsight Monitor, under *WDDM TDR Delay*
and *WDDM TDR Enabled*). Data-centre GPUs can be put in TCC mode instead.

**After a CUDA error, every call fails.**
A failed CUDA call can leave the GPU context in a bad state. Restart MATLAB
or the Python process.

**The reconstruction is blurred, doubled, shaded or full of rings.**
This is almost always the data or the geometry, not the algorithm. Compare
your image with the pictures in [](artefacts.md), and check the geometry step
by step with [](geometry-walkthrough.md).

**A Python algorithm ignores an option.**
Python prints a warning, rather than raising an error, for keyword arguments
it does not recognise. Check the spelling: for example the step size is
`lmbda`, not `lambda`.

## Using TIGRE

**Can I use TIGRE for my project, including commercial work?**
Yes. TIGRE is released under the permissive
[BSD 3-clause licence](https://github.com/CERN/TIGRE/blob/master/LICENSE.txt).
If you publish work that uses it, please [cite TIGRE](../citing.md).

**How can I contribute?**
See the [contributing guide](https://github.com/CERN/TIGRE/blob/master/contributing.md).
Bug fixes, new algorithms, scanner loaders and documentation are all
welcome.
