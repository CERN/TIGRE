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

}  // namespace tigre

#endif  // TIGRE_LAUNCH_SHAPES_HPP
