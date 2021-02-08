/*-------------------------------------------------------------------------
 *
 * MATLAB MEX  functions for TV image denoising. Check inputs and parses 
 * MATLAB data to C++ data.
 *
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
#include <CUDA/tvdenoising.hpp>
#include <CUDA/GpuIds.hpp>
#include <CUDA/gpuUtils.hpp>
/**
 * MEX gateway
 */
void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
    int maxIter;
    float lambda;
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
    if (nrhs == 0) {
        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "At least one input argumet required.");
    } else if (nrhs==1){
        maxIter=100;
        lambda=15.0f;
    } else if (nrhs==2){
        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "Only 1 TV hyperparemter inputed");
    } else if (nrhs==3 || nrhs==4){
        size_t mrows = mxGetM(prhs[1]);
        size_t ncols = mxGetN(prhs[1]);
        if (mrows!=1 || ncols !=1) {
            mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "TV parameters should be 1x1");
        }
        mrows = mxGetM(prhs[2]);
        ncols = mxGetN(prhs[2]);
        if (mrows!=1 || ncols !=1) {
            mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "TV parameters should be 1x1");
        }
        lambda= (float)(mxGetScalar(prhs[1]));
        maxIter=(int)round(mxGetScalar(prhs[2]));
    } else if (nrhs>4) {
        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "Too many imput argumets");
    }
    ////////////////////////// First input.
    // First input should be x from (Ax=b), or the image.
    mxArray const * const image = prhs[0];
    mwSize const numDims = mxGetNumberOfDimensions(image);
    
    // Image should be dim 3
    if (numDims!=3){
        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "Image is not 3D");
    }
    // Now that input is ok, parse it to C data types.
    float  *  img = static_cast<float  *>(mxGetData(image));
    // We need a float image, and, unfortunatedly, the only way of casting it is by value
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    
    //////////////
    //prepareotputs
    plhs[0] = mxCreateNumericArray(3,size_img, mxSINGLE_CLASS, mxREAL);
    float *imgout =(float*) mxGetPr(plhs[0]);
    // Allocte output image
    // call C function with the CUDA denoising
    const float spacing[3]={1,1,1};
    const long imageSize[3]={size_img[0] ,size_img[1],size_img[2] };
   
    tvdenoising(img,imgout, lambda, spacing, imageSize, maxIter, gpuids); 
    
    
    
//     memcpy(mxImgout,imgout,size_img[0] *size_img[1] *size_img[2]*sizeof(float));
    //free memory
//     free(img);
//     free(imgout);
     

}
