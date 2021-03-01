/*
/*-------------------------------------------------------------------------
 *
 * MATLAB MEX gateway for Total variation minimization via Steepest descend
 *
 * This file gets the data from MATLAB, checks it for errors and then 
 * parses it to C and calls the relevant C/CUDA fucntions.
 *
 * CODE by       Ander Biguri
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
#include <CUDA/POCS_TV.hpp>
#include <CUDA/GpuIds.hpp>
#include <CUDA/gpuUtils.hpp>
void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
///////// First check if the amount of imputs is right.    
    int maxIter;
    float alpha;
    GpuIds gpuids;
    if (nrhs==4) {
        size_t iM = mxGetM(prhs[3]);
        if (iM != 1) {
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:unknown","4th parameter must be a row vector.");
            return;
        }
        size_t uiGpuCount = mxGetN(prhs[3]);
        if (uiGpuCount == 0) {
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:unknown","4th parameter must be a row vector.");
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
    if (nrhs==1){
        maxIter=100;
        alpha=15.0f;
    } else if (nrhs==2){
       mexErrMsgIdAndTxt("minTV:mex", "Only 1 POCS hyperparemter inputed");
    } else if (nrhs==3 || nrhs==4){
     size_t mrows = mxGetM(prhs[1]);
     size_t ncols = mxGetN(prhs[1]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("minTV:mex", "POCS parameters should be 1x1");
     mrows = mxGetM(prhs[2]);
     ncols = mxGetN(prhs[2]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("minTV:mex", "POCS parameters should be 1x1");
     alpha= (float)(mxGetScalar(prhs[1]));
     maxIter=(int)floor(mxGetScalar(prhs[2])+0.5);
    } else {
       mexErrMsgIdAndTxt("minTV:mex", "Too many imput argumets");
    }
    
////////////////////////// First input.
    // First input should be x from (Ax=b), or the image.
    mxArray const * const image = prhs[0];
    mwSize const numDims = mxGetNumberOfDimensions(image);
    mwSize third_dim = 1;
    
    
    // Now that input is ok, parse it to C data types.
    float  *  img = static_cast<float  *>(mxGetData(image));
    const mwSize *size_img = mxGetDimensions(image); //get size of image    

    // Image should be dim 3
    if (numDims==3){
        third_dim = size_img[2];
    }
    
    // Allocte output image  
    const long imageSize[3]={size_img[0] ,size_img[1], third_dim };
    plhs[0] = mxCreateNumericArray(numDims, size_img, mxSINGLE_CLASS, mxREAL);
    float *imgout =(float*) mxGetPr(plhs[0]);
    
    pocs_tv(img,imgout, alpha, imageSize, maxIter, gpuids); 
}
