# Projection and backprojection

All of TIGRE is built on two GPU operators:

- `Ax` (forward projection) takes an image and simulates the projections a
  scanner would measure. It computes $A x$, where $A$ is the system matrix.
- `Atb` (backprojection) takes projections and spreads them back over the
  image. It approximates $A^T b$.

Every algorithm is a sequence of calls to these two, so understanding their
options explains most of TIGRE's behaviour.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python
```python
import numpy as np
import tigre
from tigre.utilities import sample_loader

geo = tigre.geometry_default(high_resolution=False)
angles = np.linspace(0, 2 * np.pi, 100, endpoint=False)
head = sample_loader.load_head_phantom(geo.nVoxel)

proj = tigre.Ax(head, geo, angles)                 # Siddon projector
proj_i = tigre.Ax(head, geo, angles, "interpolated")
back = tigre.Atb(proj, geo, angles)                # FDK-weighted backprojector
back_m = tigre.Atb(proj, geo, angles, "matched")
```
Both operators require `float32` arrays. The image must have shape
`geo.nVoxel` and the projections shape `(len(angles), *geo.nDetector)`.
:::

:::{tab-item} MATLAB
:sync: matlab
```matlab
geo = defaultGeometry('nVoxel', [64; 64; 64]);
angles = linspace(0, 2*pi, 100);
head = headPhantom(geo.nVoxel);

proj = Ax(head, geo, angles);                      % Siddon projector
proj_i = Ax(head, geo, angles, 'interpolated');
back = Atb(proj, geo, angles);                     % FDK-weighted backprojector
back_m = Atb(proj, geo, angles, 'matched');
```
Both operators require `single` arrays.
:::
::::

## Forward projectors

`Siddon` (the default; `ray-voxel` is an older name for it)
: Follows each ray from the source to a detector pixel and adds up the exact
  length of the ray inside each voxel it crosses, times that voxel's value.
  It is fast, and it is the projector the iterative algorithms use.

`interpolated`
: Samples the image at regular steps along each ray, using trilinear
  interpolation, and adds up the samples. The step is `geo.accuracy`, in
  voxels (default 0.5). Smaller steps are more accurate and slower; values
  between 0.1 and 0.5 are a good range. When voxels are large compared to
  detector pixels, Siddon shows discretisation artefacts that the
  interpolated projector avoids.

Demo 12
([Python](https://github.com/CERN/TIGRE/blob/master/Python/demos/d12_ProjectionOperations.py),
[MATLAB](https://github.com/CERN/TIGRE/blob/master/MATLAB/Demos/d12_ProjectionOperators.m))
compares the two and the effect of `geo.accuracy`.

## Backprojectors

`FDK` (the default)
: A voxel-driven backprojector with the distance weights used by the FDK
  algorithm. The SART family and FDK itself use it.

`matched`
: The same backprojection with weights chosen so that it is close to the
  transpose of `Ax`. The Krylov algorithms (CGLS, LSQR, LSMR and their
  variants) use it.

## The pair is not an exact adjoint

`Ax` is ray-driven and `Atb` is voxel-driven, so `Atb` is not exactly the
transpose of `Ax`, even with `matched`. Computing both operators on the fly
on the GPU, rather than storing the matrix, is what makes TIGRE fast and
lets it handle very large images, and the small mismatch is the price.

In practice:

- **For most work this does not matter.** FDK and the SART family were never
  meant to use an exact transpose, and the Krylov algorithms converge well
  with the matched backprojector for the number of iterations people use.
- **Krylov algorithms may stall after many iterations.** CGLS, LSQR and LSMR
  assume an exact transpose. If you run them for a long time, the residual
  may stop going down. Stop early, or use a method built for unmatched
  pairs.
- **AB-GMRES and BA-GMRES are built for unmatched pairs.** They do not assume
  `Atb` is the transpose of `Ax`. Their `backprojector` option set to `"FDK"`
  uses the full FDK algorithm (filtering included) as the backprojection
  step (demo 8).
- **Adjoint tests will not pass to machine precision.** If you check
  $\langle Ax, b \rangle = \langle x, A^T b \rangle$, expect the two sides to
  be close, but not equal to rounding error.

If you need an exactly matched pair, for example for a proof of convergence
or for training some learned methods, keep this mismatch in mind. To use
TIGRE's operators inside PyTorch, see the PyTorch bindings demo
([d25_Pytorch.py](https://github.com/CERN/TIGRE/blob/master/Python/demos/d25_Pytorch.py)).

## Using the operators in your own code

`Ax` and `Atb` are ordinary functions, so you can write your own algorithm on
top of them. Keep these in mind:

- The arrays must be `float32` (`single` in MATLAB). TIGRE raises an error
  rather than converting them for you.
- Each call copies data to the GPU and back. For small problems in a tight
  loop, that copy can cost more than the computation. Projecting all angles
  in one call is much faster than one angle at a time.
- Both functions accept a `gpuids` argument to choose which GPUs to use; see
  [](gpus.md).
