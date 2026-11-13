# GPUs and memory

TIGRE runs its operators on NVIDIA GPUs through CUDA. It uses every GPU it
can see by default, and it splits problems that do not fit in GPU memory into
pieces, so the size of what you can reconstruct is limited by the computer's
RAM rather than by the GPU.

## Choosing GPUs

`GpuIds` describes a set of GPUs. Create one, then pass it to any operator or
algorithm with the `gpuids` option.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python
```python
from tigre.utilities import gpu

print(gpu.getGpuNames())               # names of the installed GPUs

gpuids = gpu.GpuIds()                  # all GPUs
gpuids = gpu.GpuIds("NVIDIA RTX A6000")  # all GPUs with this name
gpuids.devices = [0, 1]                # or pick device numbers

proj = tigre.Ax(img, geo, angles, gpuids=gpuids)
img = algs.ossart(proj, geo, angles, 20, gpuids=gpuids)
```
:::

:::{tab-item} MATLAB
:sync: matlab
```matlab
getGpuNames()                          % names of the installed GPUs

gpuids = GpuIds();                     % all GPUs
gpuids = GpuIds('NVIDIA RTX A6000');   % all GPUs with this name
gpuids.devices = int32([0 1]);         % or pick device numbers

proj = Ax(img, geo, angles, 'gpuids', gpuids);
img = OS_SART(proj, geo, angles, 20, 'gpuids', gpuids);
```
:::
::::

The `CUDA_VISIBLE_DEVICES` environment variable also works, as it does for
any CUDA program: TIGRE only sees the GPUs it lists.

Demo 22 ([Python](https://github.com/CERN/TIGRE/blob/master/Python/demos/d22_ListGPUs.py),
[MATLAB](https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d22_ListGPUs.m))
shows these options.

## How the work is split

- **Forward projection** divides the angles between the GPUs. Each GPU
  computes the projections for its share of the angles.
- **Backprojection** divides the image between the GPUs, slice by slice.
- **If the image does not fit** in GPU memory, TIGRE splits it into pieces and
  processes them one after another, then puts the result together. This
  happens automatically.

The [multi-GPU paper](https://doi.org/10.1016/j.jpdc.2020.07.004) describes
the method and its performance, including reconstructions much larger than a
single GPU's memory.

Some consequences:

- **Use identical GPUs together.** TIGRE plans the split using the free
  memory of the most limited GPU and assumes all GPUs are equally fast. It
  warns you if the selected GPUs differ. With mixed GPUs, select the
  matching ones by name.
- **Keep the GPUs free.** TIGRE stops with an error if a selected GPU has less
  than half of its memory free, which usually means another program (often
  the desktop or a browser) is using it.
- **Small problems gain little from more GPUs.** Copying data to each GPU has
  a fixed cost. For small images, one GPU can be as fast as several.

## Host memory

The data you pass to TIGRE, and the arrays each algorithm keeps between
iterations, live in the computer's RAM. Krylov algorithms such as CGLS and
LSQR keep several image-sized and projection-sized arrays, so they need more
RAM than OS-SART. Hybrid LSQR keeps one image-sized and one projection-sized
vector per iteration, so its memory grows with the number of iterations.

### Pinned memory

When TIGRE splits a problem, or uses more than one GPU, it page-locks
("pins") the host arrays so that copies to and from the GPUs can run while
the GPUs compute. If you reconstruct datasets larger than your RAM and rely
on swap, or pinning fails on your system, build TIGRE without it:

- Python: set `TIGRE_NO_PINNED_MEMORY=1` when installing from source (see
  [](../install/python.md)).
- MATLAB: compile with `Compile -no_pinned_mem` (see
  [](../install/matlab.md)).

The operators are slower without pinned memory.
