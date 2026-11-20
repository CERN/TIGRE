# Tuning TIGRE for your computer

:::{warning}
This page is for advanced users. The changes below edit TIGRE's CUDA source
and need a rebuild. They give modest gains at best, and wrong values can make
TIGRE fail. Try them only if speed is critical.
:::

All the files below are in
[`Common/CUDA`](https://github.com/CERN/TIGRE/tree/master/Common/CUDA). After
each change, rebuild TIGRE (`pip install .` from the repository for Python,
`Compile` for MATLAB) and time the operation again on data the size of your
real problem.

## Skip the free-memory check (largest gain)

Before every CUDA call, TIGRE asks each GPU how much memory is free, to plan
how to split the problem ([](gpus.md)). This takes up to about 0.3 seconds
per GPU per call. For a backprojection that takes around a second, that is a
large share.

If you always run on the same, otherwise idle GPUs, you can replace the query
with a fixed number. Each `.cu` file that launches GPU work ends with a
function `checkFreeMemory`. Comment out its body and set the free memory, in
bytes, directly:

```c++
*mem_GPU_global = 11000000000;  // for example, about 11 GB
```

Use about 95% of the GPU's memory, and less if other programs (such as the
desktop) use the GPU. A value larger than the memory actually free will make
TIGRE fail. The function appears in `Siddon_projection.cu`,
`ray_interpolated_projection.cu`, `voxel_backprojection.cu`,
`voxel_backprojection2.cu`, `tv_proximal.cu`, `GD_TV.cu` and `GD_AwTV.cu`.

## Kernel launch sizes (small gains)

The constants below set how much work each GPU thread or block does. The
defaults were tuned on an NVIDIA GTX 1080 Ti. Expect at most a few per cent
from changing them, and often nothing. Change one at a time, by small steps.

| Operation | File | Constants |
|---|---|---|
| Backprojection, FDK weights | `voxel_backprojection.cu` | `PROJ_PER_KERNEL`, `VOXELS_PER_THREAD` |
| Backprojection, matched weights | `voxel_backprojection2.cu` | `PROJ_PER_KERNEL`, `VOXELS_PER_THREAD` |
| Forward projection, Siddon | `Siddon_projection.cu` | `PROJ_PER_BLOCK`, `PIXEL_SIZE_BLOCK` |
| Forward projection, interpolated | `ray_interpolated_projection.cu` | `PROJ_PER_BLOCK`, `PIXEL_SIZE_BLOCK` |
| TV denoising | `tv_proximal.cu` | `MAX_BUFFER`, `BLOCK_SIZE` |
| TV gradient steps (ASD-POCS family) | `GD_TV.cu`, `GD_AwTV.cu` | `MAX_BUFFER` |

For the forward projectors, a block has
`PIXEL_SIZE_BLOCK × PIXEL_SIZE_BLOCK × PROJ_PER_BLOCK` threads, which must not
exceed your GPU's limit (1024 on current NVIDIA GPUs). For TV denoising,
`BLOCK_SIZE`³ must not exceed it either; `MAX_BUFFER` usually has the larger
effect.

A simple way to time each operation:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python
```python
import time
import numpy as np
import tigre
from tigre.utilities import sample_loader
from tigre.utilities.im_3d_denoise import im3ddenoise

geo = tigre.geometry_default()          # use your own geometry here
angles = np.linspace(0, 2 * np.pi, 360, endpoint=False)
img = sample_loader.load_head_phantom(geo.nVoxel)
proj = tigre.Ax(img, geo, angles)

for name, op in [
    ("Ax", lambda: tigre.Ax(img, geo, angles)),
    ("Atb", lambda: tigre.Atb(proj, geo, angles)),
    ("TV", lambda: im3ddenoise(img, iter=50)),
]:
    t = time.perf_counter(); op(); print(name, time.perf_counter() - t)
```
:::

:::{tab-item} MATLAB
:sync: matlab
```matlab
geo = defaultGeometry();               % use your own geometry here
angles = linspace(0, 2*pi, 360);
img = headPhantom(geo.nVoxel);
proj = Ax(img, geo, angles);

tic; Ax(img, geo, angles); toc
tic; Atb(proj, geo, angles); toc
tic; im3DDenoise(img, 'TV', 50, 15); toc
```
:::
::::
