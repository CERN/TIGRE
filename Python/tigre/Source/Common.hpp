#ifndef _COMMON_HPP_20201017_
#define _COMMON_HPP_20201017_

#define STRINGIFY(n) #n
#define TOSTRING(n) STRINGIFY(n)
#define __HERE__ __FILE__ " (" TOSTRING(__LINE__) "): "
#define PRINT_HERE printf(__HERE__);printf
// #define PRINT_HERE (void*)0

#if defined(IS_FOR_PYTIGRE)
void mexPrintf(const char*, ...);
void mexErrMsgIdAndTxt(const char* pcTag, const char* pcMsg);
void mexWarnMsgIdAndTxt(const char* pcTag, const char* pcMsg);
#else
#include "mex.h"
#include "tmwtypes.h"
#endif  // IS_TIGRE_FOR_PYTHON
#endif  // _COMMON_HPP_20201017_
