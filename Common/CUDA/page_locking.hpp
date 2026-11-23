/*-------------------------------------------------------------------------
 *
 * Page-locking (pinning) the caller's host arrays for the copies to and from
 * the GPUs: when it is worth it (worthPageLocking, the one rule all the
 * projectors and backprojectors follow) and how (hostRegister, which leaves
 * already page-locked arrays alone and falls back to pageable memory when
 * the system will not lock that much).
 *
 *-------------------------------------------------------------------------*/
#ifndef TIGRE_PAGE_LOCKING_HPP
#define TIGRE_PAGE_LOCKING_HPP

#include <cuda_runtime_api.h>
#include <stddef.h>
#include <cstdio>
#include "TIGRE_common.hpp"

namespace tigre {

// Register the caller's host array as pinned memory, if it is not already.
// Returns true only if this call registered it, i.e. if the caller of
// hostRegister must unregister it. An array that is already page-locked is
// used as it is and left alone: registered by someone else (cudaHostRegister
// then says "already registered"), or allocated page-locked, e.g. a pinned
// PyTorch tensor (cudaHostAlloc; cudaHostRegister then says "invalid
// argument", so ask CUDA what the pointer is first).
// Any other failure (typically: the system won't lock that much memory, e.g. a volume larger than Windows lets a
// process page-lock) falls back to pageable memory with a warning: page-locking only makes the copies faster and
// asynchronous, so the operation still completes, as images larger than the GPU must.
inline bool hostRegister(void* ptr, size_t bytes) {
#if CUDART_VERSION >= 10000
    cudaPointerAttributes attributes;
    if (cudaPointerGetAttributes(&attributes, ptr) == cudaSuccess && attributes.type == cudaMemoryTypeHost)
        return false;  // already page-locked
    cudaGetLastError();  // CUDA 10 reports plain (pageable) memory as an error here
#endif
    cudaError_t err = cudaHostRegister(ptr, bytes, cudaHostRegisterPortable);
    if (err == cudaErrorHostMemoryAlreadyRegistered) {
        cudaGetLastError();  // not an error for us
        return false;
    }
    if (err != cudaSuccess) {
        cudaGetLastError();  // clear it: the caller's error check must not see it
        char msg[256];
        std::snprintf(msg, sizeof(msg), "Could not page-lock %.1f GiB of host memory (%s): using pageable memory, "
                      "slower but correct.", bytes / 1073741824.0, cudaGetErrorString(err));
        mexWarnMsgIdAndTxt("TIGRE:hostRegister", msg);
        return false;
    }
    return true;
}

// Whether to page-lock (tigre::hostRegister) the caller's image, projections or backprojection result. One rule for
// all the projectors and backprojectors: only when the image is split or more than two GPUs are used.
// Page-locking makes copies faster (24 vs 5-10 GiB/s on an RTX 4070) and asynchronous, so several GPUs can copy at
// once, but the whole array must be locked first (~16 GiB/s, slower on freshly allocated memory) and unlocked at the
// end. With one GPU and one piece every array crosses once, and the copies are either hidden behind the kernels (the
// projections, streamed in chunks; the forward projections) or a single copy that locking does not make cheaper
// overall (the image before the first kernel, the result after the last one): locking only adds time. With a split
// image the data crosses again for every piece; with several GPUs every GPU copies (the image or the projections to
// each of them, the result back from each), and pageable copies run one GPU after another. That wins from about
// three GPUs: on four Quadro RTX 6000 (512^3, 1024^3, NumPy arrays) not locking was 4-15% slower with four GPUs, but
// 5-21% faster with two, where locking took longer than the copies it sped up.
// GPU arrays are never page-locked. The forward projectors' output is never page-locked: its copies already overlap
// the kernels. The TV functions keep their own, empirical, thresholds.
inline bool worthPageLocking(int isHostRegisterSupported, bool gpuData, unsigned long splits, int deviceCount) {
    return isHostRegisterSupported && !gpuData && (splits > 1 || deviceCount > 2);
}

}  // namespace tigre

#endif  // TIGRE_PAGE_LOCKING_HPP
