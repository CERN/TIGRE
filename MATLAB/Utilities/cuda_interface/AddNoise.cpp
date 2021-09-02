/*-------------------------------------------------------------------------
 *
 * MATLAB MEX  functions for Random Number Generator. Check inputs and parses 
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
#include <CUDA/RandomNumberGenerator.hpp>
#include <CUDA/GpuIds.hpp>
#include <CUDA/gpuUtils.hpp>
/**
 * MEX gateway
 * AddNoise(Im, mu, sigma, "gpuids", gpuids);
 *   poissrnd(Im)+randn(size(Im)).*sigma + mu;
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    size_t uiLen = 0;
    float fGaussMu = 0;
    float fGaussSigma = 0;

    GpuIds gpuids;
    if (nrhs==5) {
        size_t iM = mxGetM(prhs[4]);
        if (iM != 1) {
            mexErrMsgIdAndTxt( "CBCT:MEX:RNG:unknown","5th parameter must be a row vector.");
            return;
        }
        size_t uiGpuCount = mxGetN(prhs[4]);
        if (uiGpuCount == 0) {
            mexErrMsgIdAndTxt( "CBCT:MEX:RNG:unknown","5th parameter must be a row vector.");
            return;
        }
        int* piGpuIds = (int*)mxGetData(prhs[4]);
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
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("CBCT:CUDA:RNG", "At least three input argumet required.");
    } else if (nrhs==3 || nrhs==5){
        size_t mrows = mxGetM(prhs[1]);
        size_t ncols = mxGetN(prhs[1]);
        if (mrows!=1 || ncols !=1) {
            mexErrMsgIdAndTxt("CBCT:CUDA:RNG", "2nd parameter should be 1x1");
        }
        mrows = mxGetM(prhs[2]);
        ncols = mxGetN(prhs[2]);
        if (mrows!=1 || ncols !=1) {
            mexErrMsgIdAndTxt("CBCT:CUDA:RNG", "3rd parameter should be 1x1");
        }
        fGaussMu    = (float)mxGetScalar(prhs[1]);
        fGaussSigma = (float)mxGetScalar(prhs[2]);
    } else if (nrhs>4) {
        mexErrMsgIdAndTxt("CBCT:CUDA:RNG", "Too many imput argumets");
    }
    /////////////// First input argumet.
    // First input should be an array, whose elements are lambda.
    mxArray const * const image = prhs[0];
    float* pfLambdas = static_cast<float*>(mxGetData(image));
    mwSize const numDims = mxGetNumberOfDimensions(image);  // get dim of image
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    uiLen = size_img[0];    // calculate the total length
    for (int iI = 1; iI < numDims; ++iI) {
        uiLen *= size_img[iI];
    }
    //////////////
    //prepare outputs
    // Allocte output image
    plhs[0] = mxCreateNumericArray(numDims, size_img, mxSINGLE_CLASS, mxREAL);
    float *imgout =(float*) mxGetPr(plhs[0]);
    // call CUDA rng
    poisson_gaussian_1d(pfLambdas, uiLen, fGaussMu, fGaussSigma, imgout, gpuids);
}
