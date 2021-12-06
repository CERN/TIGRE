#include <stdio.h>
#include <string.h>
#include <mex.h>
#include <CUDA/gpuUtils.hpp>

void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
	if (nrhs != 0) {
		mexErrMsgIdAndTxt("MATLAB:getGpuCount_mex", "No input requred.");
		return;
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MATLAB:getGpuCount_mex", "Too many output arguments. Returns one integer.");
		return;
	}
	int iCount = GetGpuCount();
	size_t dims[2] = {1,1};
	plhs[0] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
    *((int*)mxGetData(plhs[0])) = iCount;
}
