/*-------------------------------------------------------------------------
 *
 * Release a SYCL function's resources when it is left by an error.
 *
 * The SYCL counterpart of Common/CUDA/tigre_release_on_error.hpp (#854). In
 * the Python build an error raises (mexErrMsgIdAndTxt throws, see
 * TIGRE_common.hpp), which Cython turns into a Python RuntimeError. Without
 * this, every failed call would keep the device memory it held at that point.
 *
 * Usage, at the top of a SYCL entry point:
 *
 *     tigre::ReleaseOnError onError;
 *     onError.addQueue(&tigreSyclQueue(gpuids[dev]));   // for each device
 *     ...
 *     float* d = tigre::deviceMalloc<float>(n, gpuids[dev], tag, what);
 *     onError.add([d, q]{ sycl::free(d, *q); });
 *     ...
 *     onError.dismiss();   // just before the function's own clean-up code
 *
 * Capture by value: the actions may run after the function's locals are gone.
 * Actions run in reverse order, after waiting for the queues to finish, and
 * must never throw. On the normal path dismiss() is called and nothing runs.
 *
 * The header has the same name as the CUDA one: the SYCL build puts
 * Common/SYCL before Common/CUDA on the include path.
 *
 *-------------------------------------------------------------------------*/
#ifndef TIGRE_RELEASE_ON_ERROR_HPP
#define TIGRE_RELEASE_ON_ERROR_HPP

#include <sycl/sycl.hpp>
#include <stddef.h>
#include <functional>
#include <vector>
#include "TIGRE_common.hpp"
#include "tigre_sycl.hpp"

namespace tigre {

class ReleaseOnError {
public:
    ReleaseOnError() : active_(true) {}

    ~ReleaseOnError() {
        if (active_) release();
    }

    // Queues to wait for before releasing anything (in-flight work may still
    // use the buffers).
    void addQueue(sycl::queue* q) { queues_.push_back(q); }

    template <class F>
    void add(F action) { actions_.push_back(std::function<void()>(action)); }

    void dismiss() {
        active_ = false;
        actions_.clear();
    }

private:
    void release() {
        for (size_t i = 0; i < queues_.size(); i++) {
            try {
                queues_[i]->wait();
            } catch (...) {
            }
        }
        for (size_t i = actions_.size(); i-- > 0;) {
            try {
                actions_[i]();
            } catch (...) {
            }
        }
    }

    ReleaseOnError(const ReleaseOnError&);
    ReleaseOnError& operator=(const ReleaseOnError&);

    bool active_;
    std::vector<sycl::queue*> queues_;
    std::vector<std::function<void()> > actions_;
};

// Allocate device memory on TIGRE device iDeviceId. A failed allocation is
// reported as out of memory, with a message that starts with "out of memory"
// (see tigreNotEnoughGpuMemory), as in the CUDA build (checkMalloc), so that
// callers that retry after freeing memory recognise it.
inline void* deviceMallocBytes(size_t bytes, int iDeviceId, const char* tag, const char* what) {
    void* ptr = sycl::malloc_device(bytes, tigreSyclQueue(iDeviceId));
    if (!ptr && bytes)
        tigreNotEnoughGpuMemory(tag, what, bytes, tigreSyclFreeMemory(iDeviceId));
    return ptr;
}
template <class T>
inline T* deviceMalloc(size_t count, int iDeviceId, const char* tag, const char* what) {
    return (T*)deviceMallocBytes(count * sizeof(T), iDeviceId, tag, what);
}

}  // namespace tigre

#endif  // TIGRE_RELEASE_ON_ERROR_HPP
