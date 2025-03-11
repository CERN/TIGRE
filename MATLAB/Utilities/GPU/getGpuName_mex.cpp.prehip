#include <mex.h>
#include <CUDA/gpuUtils.hpp>

void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
    // Usage: name = getGpuName_mex(int iId)
    if (nrhs != 1) {
		mexErrMsgIdAndTxt( "MATLAB:getGpuName_mex:invalidNumInputs", "One input required.");
		return;
	} else if(nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:getGpuName_mex:maxlhs", "Too many output arguments.");
		return;
    }
    
    int iId = 0;
    if (mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt( "MATLAB:getGpuName_mex:inputNotInt", "Input must be an integer.");
        return;
    } else {
        iId = *((int*)mxGetData(prhs[0]));
    }
	int iCount = GetGpuCount();
    char* pcName = (char*)mxCalloc(128, sizeof(char));
    if (iId < iCount) {
        GetGpuName(iId, pcName);
    }
    plhs[0] = mxCreateString(pcName);
}
