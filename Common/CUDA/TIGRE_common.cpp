#if defined(IS_FOR_PYTIGRE)
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "TIGRE_common.hpp"
void mexPrintf(const char* format, ...) {
    PRINT_HERE("");
    va_list argpointer;
    va_start(argpointer, format);
    vprintf(format, argpointer);
    va_end(argpointer);
}
void mexErrMsgIdAndTxt(const char* pcTag, const char* pcMsg) {
    PRINT_HERE("%s %s\n", pcTag, pcMsg);
    exit(1);
}
void mexWarnMsgIdAndTxt(const char* pcTag, const char* pcMsg) {
    PRINT_HERE("%s %s\n", pcTag, pcMsg);
}
#endif  // IS_FOR_PYTIGRE
