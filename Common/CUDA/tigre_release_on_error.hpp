/*-------------------------------------------------------------------------
 *
 * Release a CUDA function's resources when it is left by an error.
 *
 * In the Python build a CUDA error raises (see tigreCudaError in
 * TIGRE_common.hpp), which Cython turns into a Python RuntimeError. Without
 * this, every failed call would leave behind whatever it held at that point:
 * device buffers, texture arrays, streams, pinned buffers and, worst, the
 * host registration (cudaHostRegister) of the caller's own array, which makes
 * the next call with the same array fail with "already registered".
 *
 * Usage, at the top of a CUDA entry point:
 *
 *     tigre::ReleaseOnError onError;
 *     ...acquire something...
 *     onError.add([ptr, dev]{ cudaSetDevice(dev); cudaFree(ptr); });
 *     ...
 *     onError.dismiss();   // just before the function's own clean-up code
 *
 * Capture by value: the actions may run after the function's locals are gone.
 * Actions run in reverse order, after waiting for the devices to finish, and
 * must never throw.
 *
 * On the normal path dismiss() is called and nothing runs, so the success path
 * is unchanged. In the MATLAB build errors leave through MATLAB's own
 * mexErrMsgIdAndTxt. If that jump skips destructors, nothing changes from
 * before; if it runs them, this releases what the call held, which nothing
 * else would. Either way it is safe.
 *
 *-------------------------------------------------------------------------*/
#ifndef TIGRE_RELEASE_ON_ERROR_HPP
#define TIGRE_RELEASE_ON_ERROR_HPP

#include <cuda_runtime_api.h>
#include <stddef.h>
#include <functional>
#include <vector>
#include "TIGRE_common.hpp"

namespace tigre {

class ReleaseOnError {
public:
    ReleaseOnError() : active_(true) {}

    ~ReleaseOnError() {
        if (active_) release();
    }

    // Devices to wait for before releasing anything (in-flight work may still
    // use the buffers). If none are given, the current device is used.
    void addDevice(int device) { devices_.push_back(device); }

    template <class F>
    void add(F action) { actions_.push_back(std::function<void()>(action)); }

    void dismiss() {
        active_ = false;
        actions_.clear();
    }

private:
    void release() {
        if (devices_.empty()) {
            cudaDeviceSynchronize();
        } else {
            for (size_t i = 0; i < devices_.size(); i++) {
                cudaSetDevice(devices_[i]);
                cudaDeviceSynchronize();
            }
        }
        for (size_t i = actions_.size(); i-- > 0;) {
            try {
                actions_[i]();
            } catch (...) {
            }
        }
        // Do not leave an error from the release itself to be reported by the
        // next, unrelated TIGRE call.
        cudaGetLastError();
    }

    ReleaseOnError(const ReleaseOnError&);
    ReleaseOnError& operator=(const ReleaseOnError&);

    bool active_;
    std::vector<int> devices_;
    std::vector<std::function<void()> > actions_;
};

// Register the caller's host array as pinned memory, if it is not already.
// Returns true only if this call registered it, i.e. if the caller of
// hostRegister must unregister it. An array that is already registered (for
// example pinned by the user, or by PyTorch) is used as it is and left alone.
// Any other failure is left for the caller's error check.
inline bool hostRegister(void* ptr, size_t bytes) {
    cudaError_t err = cudaHostRegister(ptr, bytes, cudaHostRegisterPortable);
    if (err == cudaErrorHostMemoryAlreadyRegistered) {
        cudaGetLastError();  // not an error for us
        return false;
    }
    return err == cudaSuccess;
}

// Check a cudaMalloc. Out of memory is reported as such, with a message that
// starts with "out of memory" (see tigreNotEnoughGpuMemory), so that it reads
// the same as TIGRE's own memory check and callers that retry after freeing
// memory recognise it. Checking each cudaMalloc matters: a failed cudaMalloc
// leaves a null pointer, and the cudaMemset that usually follows replaces the
// "out of memory" error with "invalid argument".
inline void checkMalloc(cudaError_t err, const char* tag, const char* what, size_t bytes) {
    if (err == cudaSuccess) return;
    cudaGetLastError();
    if (err == cudaErrorMemoryAllocation) {
        size_t free_bytes = 0, total_bytes = 0;
        cudaMemGetInfo(&free_bytes, &total_bytes);
        tigreNotEnoughGpuMemory(tag, what, bytes, free_bytes);
    }
    tigreCudaError(tag, "Malloc error", cudaGetErrorString(err));
}

}  // namespace tigre

#endif  // TIGRE_RELEASE_ON_ERROR_HPP
