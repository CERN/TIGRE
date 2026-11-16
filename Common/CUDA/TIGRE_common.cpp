#if defined(IS_FOR_PYTIGRE)
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <vector>
#include "TIGRE_common.hpp"
void mexPrintf(const char* format, ...) {
    PRINT_HERE("");
    va_list argpointer;
    va_start(argpointer, format);
    vprintf(format, argpointer);
    va_end(argpointer);
}
// Raised to Python as RuntimeError (Cython `except +`). This used to exit(1),
// which ended the whole Python process on any CUDA error (#853).
void mexErrMsgIdAndTxt(const char* pcTag, const char* pcMsg) {
    throw std::runtime_error(std::string(pcTag) + ": " + pcMsg);
}

// Queued and raised as Python warnings after the call, so Python shows each
// one once instead of TIGRE printing it on every call. One queue per extension
// module; calls into a module are sequential (they hold the GIL).
static std::vector<std::string>& pendingWarnings() {
    static std::vector<std::string> w;
    return w;
}
void mexWarnMsgIdAndTxt(const char* pcTag, const char* pcMsg) {
    std::string text = std::string(pcTag) + ": " + pcMsg;
    std::vector<std::string>& w = pendingWarnings();
    for (size_t i = 0; i < w.size(); i++) {
        if (w[i] == text) return;
    }
    w.push_back(text);
}
size_t tigreWarningCount() {
    return pendingWarnings().size();
}
const char* tigreWarningText(size_t i) {
    return i < pendingWarnings().size() ? pendingWarnings()[i].c_str() : "";
}
void tigreClearWarnings() {
    pendingWarnings().clear();
}
#endif  // IS_FOR_PYTIGRE
