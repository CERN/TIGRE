# Included by each CUDA extension module (.pyx). TIGRE's C++ code queues its
# warnings instead of printing them on every call (TIGRE_common.cpp); this
# raises them as Python warnings, so Python shows each one once and they can be
# filtered like any other warning. Each module has its own queue.
import warnings as _tigre_py_warnings

cdef extern from "TIGRE_common.hpp":
    size_t tigreWarningCount()
    const char* tigreWarningText(size_t i)
    void tigreClearWarnings()


cdef _raise_tigre_warnings():
    cdef size_t i
    cdef size_t n = tigreWarningCount()
    messages = [(<bytes>tigreWarningText(i)).decode("utf-8", "replace") for i in range(n)]
    tigreClearWarnings()
    for message in messages:
        # stacklevel=2: the caller of tigre.Ax/Atb/..., not TIGRE's wrapper.
        _tigre_py_warnings.warn(message, RuntimeWarning, stacklevel=2)
