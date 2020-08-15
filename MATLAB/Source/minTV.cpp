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





#include "tmwtypes.h"
#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "POCS_TV.hpp"
#include <string.h>
// #include <time.h>
void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
///////// First check if the amount of imputs is rigth.    
    int maxIter;
    float alpha;
    if (nrhs==1){
        maxIter=100;
        alpha=15.0f;
    }
    if (nrhs==2){
       mexErrMsgIdAndTxt("minTV:mex", "Only 1 POCS hyperparemter inputed");
    }
    if (nrhs>3){
       mexErrMsgIdAndTxt("minTV:mex", "Too many imput argumets");
    }
    if (nrhs==3){
     size_t mrows = mxGetM(prhs[1]);
     size_t ncols = mxGetN(prhs[1]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("minTV:mex", "POCS parameters shoudl be 1x1");
     mrows = mxGetM(prhs[2]);
     ncols = mxGetN(prhs[2]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("minTV:mex", "POCS parameters shoudl be 1x1");
     alpha= (float)(mxGetScalar(prhs[1]));
     maxIter=(int)floor(mxGetScalar(prhs[2])+0.5);
    }
    
////////////////////////// First input.
    // First input should be x from (Ax=b), or the image.
    mxArray const * const image = prhs[0];
    mwSize const numDims = mxGetNumberOfDimensions(image);
    
    // Image should be dim 3
    if (numDims!=3){
      //  mexErrMsgIdAndTxt("minTV:mex", "Image is not 3D");
    }
    // Now that input is ok, parse it to C data types.
    float  *  img = static_cast<float  *>(mxGetData(image));
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    

    
    // Allocte output image  
    const long imageSize[3]={size_img[0] ,size_img[1],size_img[2] };
    plhs[0] = mxCreateNumericArray(3,size_img, mxSINGLE_CLASS, mxREAL);
    float *imgout =(float*) mxGetPr(plhs[0]);
    
    
    pocs_tv(img,imgout, alpha, imageSize, maxIter); 
    

    
}