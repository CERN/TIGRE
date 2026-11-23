/*-------------------------------------------------------------------------
 *
 * Launch shapes of the projectors and backprojectors: how a problem is laid
 * out in thread blocks and launches. They depend on the problem's shape only
 * (detector size, slices), never on its values, so any shape gives the same
 * result; they exist because the shapes that suit 3D problems leave most of
 * the GPU idle on thin ones (fan beam: one detector row, one slice).
 *
 *-------------------------------------------------------------------------*/
#ifndef TIGRE_LAUNCH_SHAPES_HPP
#define TIGRE_LAUNCH_SHAPES_HPP

#include <math.h>
#include <stddef.h>

namespace tigre {

// Thread-block shape (u, v) of the forward projectors, whose blocks are (u, v, angles): base x base detector
// pixels, except on a thin detector. With one row (fan beam) a base x base tile left base-1 of every base threads
// with no pixel, so the threads go along the long side instead, keeping about base*base of them per angle.
// The kernels compute (u, v) from blockIdx/blockDim and check the bounds, so any shape gives the same result.
inline void detectorBlockShape(unsigned long nU, unsigned long nV, int base, int* divU, int* divV) {
    const int n = base * base;
    if (nV < (unsigned long)base) {
        *divV = nV > 0 ? (int)nV : 1;
        *divU = (n + *divV - 1) / *divV;
    } else if (nU < (unsigned long)base) {
        *divU = nU > 0 ? (int)nU : 1;
        *divV = (n + *divU - 1) / *divU;
    } else {
        *divU = base;
        *divV = base;
    }
}

// Angles per forward-projection launch: perBlock (one block of angles), or on a thin detector (fewer rows than
// perBlock, e.g. fan beam) up to maxPerLaunch: there one block of angles is too little work to keep the GPU busy and
// every launch is followed by a wait on the host. The two projection buffers stay within ~4 MB each. Measured on an
// RTX 4070: fan-beam Ax 2-4x faster; on 3D detectors batching did not help, so they keep one block per launch.
inline int projectionsPerLaunch(unsigned long nV, size_t bytesPerProjection, int perBlock, int maxPerLaunch) {
    if (nV >= (unsigned long)perBlock) return perBlock;
    const size_t budget = (size_t)4 << 20;
    size_t k = budget / ((size_t)perBlock * (bytesPerProjection > 0 ? bytesPerProjection : 1));
    const size_t kmax = (size_t)(maxPerLaunch / perBlock);
    if (k < 1) k = 1;
    if (k > kmax) k = kmax;
    return (int)k * perBlock;
}

// Voxels per thread along z for the cone backprojectors, which give each thread a column of voxels and, for speed,
// compute the whole column even past the end of the image. On a thin image (fan beam: one slice) most of a long
// column would be computed for nothing, so use the longest power of two, up to maxLen, that fits in the image.
inline int voxelColumnLength(unsigned long nVoxelZ, int maxLen) {
    int len = 1;
    while (len * 2 <= maxLen && (unsigned long)(len * 2) <= nVoxelZ) len *= 2;
    return len;
}

// Whether the cone backprojectors upload all the projections into one texture at once, instead of in chunks of
// PROJ_PER_KERNEL (32), each copied while the GPU works on the previous one. On a thin detector (fewer rows than
// thinRows: fan beam has one) a chunk is a few tens of KB and copies in microseconds, so there is nothing to hide,
// while every chunk costs driver calls (the copy, re-creating the texture, waits): ~0.1 ms each on an RTX 4070, most
// of a small fan-beam backprojection. Up to maxProjections, within the deepest 3D texture every CUDA GPU supports
// (2048). The caller also checks that they fit in GPU memory next to the image.
inline bool projectionsInOneUpload(unsigned long nDetecV, unsigned long nalpha, unsigned long thinRows = 9, unsigned long maxProjections = 2048) {
    return nDetecV < thinRows && nalpha <= maxProjections;
}

// Steps per pass of the Siddon forward projectors' walk along a ray (1, 2 or 4: the kernels are compiled for each,
// and give the same sums for any). raySpacing: how far apart neighbouring rays are where they cross the image, in
// voxels; voxels: the size of the image, or of the piece of it in one texture. More steps per pass put more texture
// fetches in flight, which hides their latency as long as the voxels each ray reads next are still in the cache. They
// are when neighbouring rays are at most about a voxel apart (the usual detector, a pixel per voxel): four steps.
// When the rays are further apart, the extra fetches in flight evict them first (Nsight Compute, 1024^3 image, rays 4
// voxels apart: L1 hit rate 57% with one step per pass, 41% with four), and four steps were up to 1.5x slower than
// one: two steps, and one when the image is also large. Measured on an RTX 4070 and a Quadro RTX 6000, cone beam
// 256^3..1024^3 and parallel beam 512^3, 1024^3, rays 0.5..6 voxels apart.
inline int siddonStepsPerPass(double raySpacing, size_t voxels) {
    if (raySpacing < 1.25) return 4;
    if (raySpacing * cbrt((double)voxels) < 2700.0) return 2;
    return 1;
}

}  // namespace tigre

#endif  // TIGRE_LAUNCH_SHAPES_HPP
