
#ifndef TVDENOISE
#define TVDENOISE
#include "mex.h"
#include "tmwtypes.h"
void tvdenoising(const float* src, float* dst, float lambda,
                 const float* spacing,const long* image_size, int maxIter);

#endif