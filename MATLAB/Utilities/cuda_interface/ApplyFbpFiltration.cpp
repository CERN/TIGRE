/*-------------------------------------------------------------------------
 *
 * MATLAB MEX  functions for convolution. Check inputs and parses 
 * MATLAB data to C++ data.
 *
 *
 * CODE by       Tomoyuki SADAKANE
 *
---------------------------------------------------------------------------
---------------------------------------------------------------------------
Copyright (c) 2015, University of Bath and CERN- European Organization for 
Nuclear Research
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
 ---------------------------------------------------------------------------

Contact: tigre.toolbox@gmail.com
Codes  : https://github.com/CERN/TIGRE
--------------------------------------------------------------------------- 
 */

#include <math.h>
#include <string.h>
#include <tmwtypes.h>
#include <mex.h>
#include <matrix.h>
#include <CUDA/FbpFiltration.hpp>
#include <CUDA/GpuIds.hpp>
#include <CUDA/gpuUtils.hpp>
/**
 * MEX gateway
 * ApplyFbpFiltration(mat2dProj, mat1dFlt, "gpuids", gpuids);
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    GpuIds gpuids;
    if (nrhs==4) {
        size_t iM = mxGetM(prhs[3]);
        if (iM != 1) {
            mexErrMsgIdAndTxt( "TIGRE:MEX:ApplyFbpFiltration:unknown","5th parameter must be a row vector.");
            return;
        }
        size_t uiGpuCount = mxGetN(prhs[3]);
        if (uiGpuCount == 0) {
            mexErrMsgIdAndTxt( "TIGRE:MEX:ApplyFbpFiltration:unknown","5th parameter must be a row vector.");
            return;
        }
        int* piGpuIds = (int*)mxGetData(prhs[3]);
        gpuids.SetIds(uiGpuCount, piGpuIds);
    } else {
        int iGpuCount = GetGpuCount();
        int* piDev = (int*)malloc(iGpuCount * sizeof(int));
        for (int iI = 0; iI < iGpuCount; ++iI) {
            piDev[iI] = iI;
        }
        gpuids.SetIds(iGpuCount, piDev);
        free(piDev); piDev = 0;
    }
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("TIGRE:MEX:ApplyFbpFiltration", "At least two input argumet required.");
    } else if (nrhs==2 || nrhs==4){
        size_t mrows = mxGetM(prhs[0]);
        size_t ncols = mxGetN(prhs[0]);
        if (mrows ==1 && ncols ==1) {
            mexErrMsgIdAndTxt("TIGRE:MEX:ApplyFbpFiltration", "1st parameter should not be 1x1");
        }
        mrows = mxGetM(prhs[1]);
        ncols = mxGetN(prhs[1]);
        if (mrows==1 || ncols !=1) {
            mexErrMsgIdAndTxt("TIGRE:MEX:ApplyFbpFiltration", "2nd parameter should be Mx1 (M>1)");
        }
    } else if (nrhs>4) {
        mexErrMsgIdAndTxt("TIGRE:MEX:ApplyFbpFiltration", "Too many imput argumets");
    }
    /////////////// First input argumet. Projection
    mxArray const * const proj = prhs[0];
    float* pfProj = static_cast<float*>(mxGetData(proj));
    const mwSize numDimsProj = mxGetNumberOfDimensions(proj);  // get dim of proj
    const mwSize *size_proj= mxGetDimensions(proj); //get size of proj
    // printf("numDimsProj = %d\n", numDimsProj);
    // for (int iI = 0; iI < numDimsProj; ++iI) {
    //    printf("size_proj[%d] = %d\n", iI, size_proj[iI]);
    // }
    /////////////// Second input argumet. Filter
    const mxArray* filter = prhs[1];
    float* pfFilter = static_cast<float*>(mxGetData(filter));
    const mwSize numDimsFilter = mxGetNumberOfDimensions(filter);  // get dim of filter
    const mwSize *size_filter= mxGetDimensions(filter); //get size of filter
    // printf("numDimsFilter = %d\n", numDimsFilter);
    // for (int iI = 0; iI < numDimsFilter; ++iI) {
    //    printf("size_filter[%d] = %d\n", iI, size_filter[iI]);
    // }
    if (size_filter[0] != size_proj[0]) {
        mexErrMsgIdAndTxt("TIGRE:MEX:ApplyFbpFiltration", "Width of projection must be equal to filter.");
    }
    //////////////
    //prepare outputs
    // Allocate output projection
    plhs[0] = mxCreateNumericArray(numDimsProj, size_proj, mxSINGLE_CLASS, mxREAL);
    float *pfProjOut =(float*) mxGetPr(plhs[0]);
    // call CUDA filtering
    apply_filtration(pfProj, size_proj[0], size_proj[1], pfFilter, pfProjOut, gpuids);
}
