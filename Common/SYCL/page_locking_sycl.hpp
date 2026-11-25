/*-------------------------------------------------------------------------
 *
 * Page-locking the caller's host arrays for the copies to and from the
 * devices, as Common/CUDA/page_locking.hpp does in the CUDA build: when it is
 * worth it (worthPageLocking, the same rule) and how (hostRegister, which
 * leaves arrays that are already device-accessible alone and falls back to
 * pageable memory when the system will not lock that much).
 *
 * SYCL: sycl::ext::oneapi::experimental::prepare_for_device_copy (extension
 * sycl_ext_oneapi_copy_optimize) for every context of the devices, released
 * with release_from_device_copy.
 *
 * The SYCL counterpart of Common/CUDA/page_locking.hpp.
 *
 *-------------------------------------------------------------------------*/
#ifndef TIGRE_PAGE_LOCKING_SYCL_HPP
#define TIGRE_PAGE_LOCKING_SYCL_HPP

#include <sycl/sycl.hpp>
#include <stddef.h>
#include <cstdio>
#include <vector>
#include "TIGRE_common.hpp"
#include "tigre_sycl.hpp"

namespace tigre {

// The SYCL contexts of TIGRE devices devs[0..ndev-1], each once (usually one: DPC++
// gives the devices of a platform one default context).
inline std::vector<sycl::context> deviceContexts(const int* devs, int ndev) {
    std::vector<sycl::context> contexts;
    for (int d = 0; d < ndev; d++) {
        sycl::context c = tigreSyclQueue(devs[d]).get_context();
        bool seen = false;
        for (size_t i = 0; i < contexts.size(); i++) seen = seen || contexts[i] == c;
        if (!seen) contexts.push_back(c);
    }
    return contexts;
}

// Whether host arrays can be page-locked for these devices (cudaDevAttrHostRegisterSupported in the
// CUDA build): GPUs only, with sycl_ext_oneapi_copy_optimize. A CPU device reads host memory as it is.
inline int hostRegisterSupported(const int* devs, int ndev) {
#ifdef SYCL_EXT_ONEAPI_COPY_OPTIMIZE
    for (int d = 0; d < ndev; d++)
        if (!tigreSyclDevices()[devs[d]].is_gpu()) return 0;
    return ndev > 0;
#else
    return 0;
#endif
}

// Page-lock the caller's host array for the devices' copies, if it is not already device-accessible.
// Returns true only if this call locked it, i.e. if the caller must release it (hostUnregister).
// An array that SYCL already knows (a USM allocation) is used as it is and left alone.
// A failure (typically: the system won't lock that much memory) falls back to pageable memory
// with a warning: page-locking only makes the copies faster, so the operation still completes.
inline bool hostRegister(void* ptr, size_t bytes, const int* devs, int ndev) {
#ifdef SYCL_EXT_ONEAPI_COPY_OPTIMIZE
    std::vector<sycl::context> contexts = deviceContexts(devs, ndev);
    for (size_t i = 0; i < contexts.size(); i++)
        if (sycl::get_pointer_type(ptr, contexts[i]) != sycl::usm::alloc::unknown)
            return false;  // already device-accessible
    size_t done = 0;
    try {
        for (; done < contexts.size(); done++)
            sycl::ext::oneapi::experimental::prepare_for_device_copy(ptr, bytes, contexts[done]);
    } catch (const sycl::exception& e) {
        for (size_t i = 0; i < done; i++) {
            try { sycl::ext::oneapi::experimental::release_from_device_copy(ptr, contexts[i]); } catch (...) {}
        }
        char msg[256];
        std::snprintf(msg, sizeof(msg), "Could not page-lock %.1f GiB of host memory (%s): using pageable memory, "
                      "slower but correct.", bytes / 1073741824.0, e.what());
        mexWarnMsgIdAndTxt("TIGRE:hostRegister", msg);
        return false;
    }
    return true;
#else
    (void)ptr; (void)bytes; (void)devs; (void)ndev;
    return false;
#endif
}

// Release an array that hostRegister locked. Never throws.
inline void hostUnregister(void* ptr, const int* devs, int ndev) {
#ifdef SYCL_EXT_ONEAPI_COPY_OPTIMIZE
    std::vector<sycl::context> contexts = deviceContexts(devs, ndev);
    for (size_t i = 0; i < contexts.size(); i++) {
        try { sycl::ext::oneapi::experimental::release_from_device_copy(ptr, contexts[i]); } catch (...) {}
    }
#else
    (void)ptr; (void)devs; (void)ndev;
#endif
}

// Whether to page-lock (tigre::hostRegister) the caller's image, projections or backprojection result: the rule of
// the CUDA build (Common/CUDA/page_locking.hpp), only when the image is split or more than two devices are used.
inline bool worthPageLocking(int isHostRegisterSupported, bool gpuData, unsigned long splits, int deviceCount) {
    return isHostRegisterSupported && !gpuData && (splits > 1 || deviceCount > 2);
}

}  // namespace tigre

#endif  // TIGRE_PAGE_LOCKING_SYCL_HPP
