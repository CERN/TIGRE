#ifndef _COMMON_HPP_20201017_
#define _COMMON_HPP_20201017_

#define STRINGIFY(n) #n
#define TOSTRING(n) STRINGIFY(n)
#define __HERE__ __FILE__ " (" TOSTRING(__LINE__) "): "
#define PRINT_HERE printf(__HERE__);printf
// #define PRINT_HERE (void*)0

#if defined(IS_FOR_PYTIGRE)
#ifndef IS_FOR_MATLAB_TIGRE
    #define IS_FOR_MATLAB_TIGRE 0
#endif  // IS_FOR_MATLAB_TIGRE
#include <stddef.h>
#include <stdexcept>
#include <string>
void mexPrintf(const char*, ...);
// Python build: errors throw std::runtime_error, which Cython's `except +`
// turns into a Python RuntimeError (#853). Warnings are queued and raised as
// Python warnings by the Cython layer after each call (tigreWarningCount etc.).
void mexErrMsgIdAndTxt(const char* pcTag, const char* pcMsg);
void mexWarnMsgIdAndTxt(const char* pcTag, const char* pcMsg);
size_t tigreWarningCount();
const char* tigreWarningText(size_t i);
void tigreClearWarnings();
// Report a failed CUDA call. `where` is what TIGRE was doing, `what` the CUDA
// error string.
inline void tigreCudaError(const char* tag, const char* where, const char* what) {
    throw std::runtime_error(std::string(tag) + ": " + where + ": " + what);
}
// Code that only runs in the MATLAB build. Used for the cudaDeviceReset() some
// error macros do before raising: in Python, where the error can now be caught,
// it would destroy every other GPU allocation in the process (e.g. PyTorch's).
#define TIGRE_MATLAB_ONLY(code)
#else
#ifndef IS_FOR_MATLAB_TIGRE
    #define IS_FOR_MATLAB_TIGRE 1
#endif  // IS_FOR_MATLAB_TIGRE
#include "mex.h"
#include "tmwtypes.h"
// MATLAB build: unchanged behaviour, print where it failed and raise a MATLAB
// error.
inline void tigreCudaError(const char* tag, const char* where, const char* what) {
    mexPrintf("%s \n", where);
    mexErrMsgIdAndTxt(tag, what);
}
#define TIGRE_MATLAB_ONLY(code) code
#endif  // IS_TIGRE_FOR_PYTHON

#include <stdio.h>
// Not even the smallest piece of the problem fits in the free GPU memory:
// raise a clear error before allocating anything. Before, the split
// calculation subtracted past zero (size_t) and Atb then divided by zero,
// killing the process. The message starts with "out of memory", so callers
// that retry after freeing memory (TIGRE's PyTorch bindings) recognise it.
inline void tigreNotEnoughGpuMemory(const char* tag, const char* what, size_t needed, size_t available) {
    char msg[512];
    snprintf(msg, sizeof(msg),
             "out of memory: %s needs at least %llu MB of free GPU memory on each GPU, "
             "but only about %llu MB are free. Free some GPU memory (in PyTorch: "
             "torch.cuda.empty_cache()) or use a smaller problem.\n",
             what, (unsigned long long)(needed >> 20) + 1, (unsigned long long)(available >> 20));
    mexErrMsgIdAndTxt(tag, msg);
}
#endif  // _COMMON_HPP_20201017_
