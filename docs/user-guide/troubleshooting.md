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

**`texture memory fail` or `invalid device function`.**
TIGRE was not built for your GPU. Rebuild it on the computer you use, with a
CUDA version that supports your GPU.

**`One (or more) of your GPUs is being heavily used by another program`.**
TIGRE needs at least half of each selected GPU's memory to be free. Close
programs using the GPU, or select another GPU; see [](gpus.md).

**`the launch timed out and was terminated` or `unspecified launch failure`
on large images.**
The operating system stopped a GPU kernel because it ran too long, which
happens on GPUs that also drive a display. On Windows, raise the TDR delay or
turn TDR off (for example in NVIDIA Nsight Monitor, under *WDDM TDR Delay*
and *WDDM TDR Enabled*). Data-centre GPUs can be put in TCC mode instead.

**After a CUDA error, every call fails.**
A failed CUDA call can leave the GPU context in a bad state. Restart MATLAB
or the Python process.

**The reconstruction is blurred, doubled or shaded.**
This is almost always the geometry, not the algorithm. Work through the
checklist at the end of [](geometry.md).

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
