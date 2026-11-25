---
og:description: "Use TIGRE with other tools: PyTorch operators with gradients for learned reconstruction, CIL and ASTRA."
myst:
  html_meta:
    description: "Use TIGRE with other tools: PyTorch operators with gradients for learned reconstruction, CIL and ASTRA."
---

# Using TIGRE with other tools

TIGRE works on plain arrays (NumPy in Python, `single` arrays in MATLAB), so
it fits into most imaging pipelines. This page covers the tools people most
often combine it with.

## Arrays in and out

Everything TIGRE takes or returns is an ordinary array:

- **Python:** `float32` NumPy arrays. Images have shape `geo.nVoxel`, indexed
  `img[z, y, x]`; projections have shape `(len(angles), *geo.nDetector)`,
  indexed `proj[angle, v, u]`.
- **MATLAB:** `single` arrays, indexed `img(x, y, z)` and `proj(v, u, angle)`.

So any reader that produces such an array (`tifffile`, `h5py`, `imageio`,
`pydicom`, MATLAB's `imread` or `h5read`) can feed TIGRE, and the results
can be saved with the same tools. Remember that TIGRE expects line integrals,
not raw intensities: see [](loading-data.md). For scanners that TIGRE reads
directly (Nikon, Bruker, YXLON, Diondo, Varian, DXchange, DICOM), the loaders
on that page also build the geometry for you.

## PyTorch

TIGRE's projectors are available as PyTorch modules with gradients, so they
can sit inside a neural network, for example in unrolled or learned iterative
reconstruction. Demo 25 trains a small learned-gradient network this way.

```python
import numpy as np
import tigre
from tigre.utilities.pytorch_bindings import create_pytorch_operator

geo = tigre.geometry(mode="fan")
angles = np.linspace(0, np.pi, 200)
A, At = create_pytorch_operator(geo, angles)   # torch.nn.Module objects

sino = A(volume)    # volume: tensor whose last three dimensions are geo.nVoxel
back = At(sino)     # sino: tensor whose last three dimensions are (angles, v, u)
```

- Any number of leading dimensions (batch, channels) is allowed; each one is
  projected separately.
- CUDA tensors stay on their GPU: TIGRE reads them and writes its result
  there, with no copy through host memory (see [GPU arrays](#gpu-arrays)).
  Tensors on the CPU go through NumPy, on all GPUs.
- `create_pytorch_operator(geo, angles, gpuids)` chooses which GPUs TIGRE
  uses. For CUDA tensors, leave it out: naming other GPUs, or several, makes
  TIGRE copy the data through host memory to use them.
- Gradients flow through both modules: the gradient of `A` is computed with
  TIGRE's `matched` backprojector and the gradient of `At` with the forward
  projector. As TIGRE's projector pairs are not exact transposes (see
  [](projectors.md)), these gradients are close to, but not exactly, the true
  ones. `create_pytorch_operator(..., backprojection_type="FDK")` uses FDK
  weights instead, as earlier versions did.
- If TIGRE runs out of GPU memory, the bindings give the memory PyTorch keeps
  in its cache back to CUDA (`torch.cuda.empty_cache()`) and try once more.

The bindings are in Python only and need `torch` installed.

## GPU arrays

`tigre.Ax` and `tigre.Atb` also take a PyTorch CUDA tensor or a
[CuPy](https://cupy.dev) array, float32, and return the same kind of array on
the same GPU. The data does not go through host memory, so this is faster
when the data is on the GPU anyway, for example between the layers of a
network.

```python
x = torch.from_numpy(img).cuda()
proj = tigre.Ax(x, geo, angles)     # a CUDA tensor, on x's GPU
back = tigre.Atb(proj, geo, angles)
```

- This uses one GPU: the one the array is on. If `gpuids` names other GPUs, or
  several, the array is copied to the host, the usual multi-GPU code runs, and
  the result is copied back; the result is the same.
- The image or the projections must fit on that GPU, with room for a copy of
  the input (TIGRE reads it through a CUDA texture). The output is a new array.
- TIGRE waits for the work that wrote the input on PyTorch's or CuPy's current
  stream, and has finished writing the result when the call returns.
- The algorithms (`tigre.algorithms`) still take NumPy arrays only.

## CIL (Core Imaging Library)

[CIL](https://tomographicimaging.github.io/CIL/) can use TIGRE as its GPU
backend. Its `cil.plugins.tigre` plugin provides a `ProjectionOperator` and
an `FBP` processor built on TIGRE's projectors, and `cil.recon.FDK` uses TIGRE
for the backprojection. CIL converts its own geometry description to TIGRE's,
so you describe the scan once in CIL and can then use CIL's optimisation
framework (for example PDHG or FISTA with any of its regularisers) with
TIGRE's projectors underneath. See the
[CIL plugin documentation](https://tomographicimaging.github.io/CIL/nightly/plugins).

TIGRE needs to be installed in the same environment as CIL for the plugin to
work.

## ASTRA Toolbox

There is no automatic converter between TIGRE and
[ASTRA](https://astra-toolbox.com/) geometries, but the two describe the same
scans, so translating by hand is straightforward. Things to watch:

| | TIGRE | ASTRA |
|---|---|---|
| Units | Millimetres everywhere (`DSO`, `DSD`, `dVoxel`, `dDetector`) | Detector spacing and distances are in units of the volume's voxel size (1 by default) |
| Source-to-origin distance | `DSO` | `source_origin` |
| Origin-to-detector distance | `DSD - DSO` | `origin_det` |
| Source position at angle 0 | On the +x axis | On the −y axis (ASTRA's `cone` geometry places the source at `(sin θ, −cos θ, 0)` times the distance) |
| Arbitrary geometries | Per-projection `offOrigin`, `offDetector`, `rotDetector`, `DSO`, `DSD` and three Euler angles per projection | `cone_vec` / `parallel3d_vec`: source, detector centre and two detector pixel vectors per projection |

Because the axis and angle conventions differ, check any translation by
projecting the same simple phantom (for example
`tigre.utilities.sample_loader.load_head_phantom`) with both tools and
comparing the projections before reconstructing real data.

## MATLAB and Python together

The MATLAB and Python versions of TIGRE use the same CUDA code and the same
geometry, but order the axes differently (see
[](geometry.md#python-and-matlab-order-the-axes-differently)). To move data
between them, save it in a format both read (for example HDF5 or `.mat` via
`scipy.io`) and reorder the axes: a MATLAB image `img(x, y, z)` becomes
`img.transpose(2, 1, 0)` in Python (`[z, y, x]`), and MATLAB projections
`proj(v, u, angle)` become `proj.transpose(2, 0, 1)` (`[angle, v, u]`).
