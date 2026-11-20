/*-------------------------------------------------------------------------
 *
 * Gradient descent on a TV-type functional, on one or several GPUs, for images
 * that may not fit on them: the GPU bookkeeping shared by minTV (GD_TV.cu) and
 * PICCS (PICCS.cu).
 *
 * The functional itself is not here: each caller gives a function that computes
 * the step direction on one GPU, for its piece of the image. This code splits the
 * image in z between GPUs (and in more pieces if it does not fit), keeps halo
 * slices around each piece and swaps them between GPUs, normalises the direction
 * by its norm over the whole image (summed in double), and steps:
 *     x <- x - alpha * direction / ||direction||
 * Extra input volumes (e.g. PICCS's prior) are split and haloed like the image.
 *
 * Compiled by nvcc only (it uses CUDA types); the MEX/Cython front ends include
 * GD_TV.hpp / PICCS.hpp.
 * ---------------------------------------------------------------------------
 * This file is part of the TIGRE Toolbox. License: BSD, see LICENSE.
 * Contact: tigre.toolbox@gmail.com   Codes: https://github.com/CERN/TIGRE
 * ---------------------------------------------------------------------------
 */
#ifndef TV_GRADIENT_DESCENT_HPP
#define TV_GRADIENT_DESCENT_HPP
#include "cuda_to_hip.h"
#include "GpuIds.hpp"

// One GPU's piece of the image, as the step function sees it. Every buffer has the same layout: `depth` slots
// (z) of rows*cols voxels. The image's slices are the slots [zlo,zhi): beyond them are halo slots holding the
// neighbouring piece's slices, or nothing at the first and the last slice of the whole image (treat voxels outside
// [zlo,zhi) like voxels outside the image in x and y). The piece itself, which is what the step updates, starts
// at element interior_offset and has interior_pixels voxels.
struct TVPiece {
    const float* image;            // the current image
    const float* const* extra;     // the extra input volumes, in the order given to tv_gradient_descent
    float* direction;              // output: the step direction (only the interior is used)
    long depth, zlo, zhi, rows, cols;
    size_t interior_offset, interior_pixels;
    cudaStream_t stream;           // launch on this stream
};

// Computes piece.direction from piece.image (and piece.extra). params: the caller's own parameters.
typedef void (*TVStepDirection)(const TVPiece& piece, const void* params);

// extra: nExtra volumes of the image's size (host memory). what: the name in out-of-memory messages.
void tv_gradient_descent(const float* img, const float* const* extra, int nExtra,
                         TVStepDirection step, const void* params,
                         float* dst, float alpha, const long* image_size, int maxIter, const GpuIds& gpuids,
                         const char* what);

#endif
