# Siddon ray tracing on GPUs

*Summary of the technical note "Siddon–Jacobs ray tracing on GPUs: a technical
note on a consistent and accurate implementation", Ander Biguri, November 2026,
revised for TIGRE 3.3.1.*
**{download}`Download the full note (PDF) <siddon_note.pdf>`**

## What the projector does

TIGRE's default forward projector, `Siddon`, computes each projection value
as the exact length of the ray inside every voxel it crosses, times that
voxel's value. Siddon's method works along the ray: the ray crosses three
families of equally spaced planes (the voxel boundaries in x, y and z), so
the lengths inside voxels are the distances between consecutive crossings.
Jacobs et al. made it faster with an incremental walk: start in the first
voxel, then repeatedly step to whichever neighbour the ray reaches first,
adding a constant to the crossing parameter of that direction. On a GPU every
ray is one thread, and the image is read through a texture.

## Three issues with the old kernel

Re-examining the kernel for issue
[#859](https://github.com/CERN/TIGRE/issues/859) (wrong projections when
the image is split into slabs) found three problems with independent causes:

| Issue | Cause | Change |
|---|---|---|
| **Consistency** | The first voxel was taken from a midpoint, and the crossings from plane limits. For an exactly characterised set of rays, those that enter the image on a voxel edge, the two disagree and the ray is integrated one voxel off. This was the cause of #859, and it also affects whole images. | The first voxel comes from the entry point and the ray direction, and the crossings from the current voxel (the initialisation of Amanatides and Woo). |
| **Accuracy** | Each crossing was found by adding a step to the previous one, so rounding errors in single precision grew along the ray. | Each crossing is computed directly from its plane, as in Siddon's original method, with one multiply-add. |
| **Speed** | The texture read was inside the branches of the step, there were conversions between integers and floats, and indices were 64-bit on Linux. | One read before the branch, voxel centres kept as floats, 32-bit indices. |

The causes are unrelated, and the speed could have been gained while keeping
Jacobs' walk. But the changes fit together: keeping the voxel as a
floating-point centre, done for speed, is what makes the direct crossings
free, and deriving the crossings from the current voxel is what makes the
start consistent. So one revised walk addresses all three at no extra cost.

## Results in TIGRE

- **Consistent for every ray**, including rays that enter on voxel edges and
  images split across slabs or GPUs.
- **About 100 times more accurate** projections, compared with exact line
  integrals.
- **1.4–6.5 times faster**, depending on the GPU and platform; most of the
  gain comes from how the loop is arranged on the GPU, not from the algorithm.
- **TIGRE 3.3.1** goes further: the walk takes up to four steps per pass,
  without branches, so that several texture reads are in flight at once. This
  makes the kernel another 1.1–1.2 times faster, with the same results. The
  number of steps is chosen from the ray spacing and the image size
  ([#913](https://github.com/CERN/TIGRE/pull/913),
  [#916](https://github.com/CERN/TIGRE/pull/916)).

The full note derives which rays the old start got wrong, analyses the
rounding errors, gives the algorithm side by side with Jacobs', separates
the effect of each change with an ablation on several GPUs, and compares the
projector with other implementations.

## Where it is in the code

The cone-beam kernel is `Common/CUDA/Siddon_projection.cu`, the parallel-beam
one `Common/CUDA/Siddon_projection_parallel.cu`. The changes were made in
[#865](https://github.com/CERN/TIGRE/pull/865) and revised for 3.3.1 in #913
and #916. To use the projector, see [](../../user-guide/projectors.md); it is
the default for `Ax`.
