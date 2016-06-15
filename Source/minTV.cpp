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
       mexErrMsgIdAndTxt("err", "Only 1 POCS hyperparemter inputed");
    }
    if (nrhs>3){
       mexErrMsgIdAndTxt("err", "Too many imput argumets");
    }
    if (nrhs==3){
     size_t mrows = mxGetM(prhs[1]);
     size_t ncols = mxGetN(prhs[1]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("err", "POCS parameters shoudl be 1x1");
     mrows = mxGetM(prhs[2]);
     ncols = mxGetN(prhs[2]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("err", "POCS parameters shoudl be 1x1");
     alpha= (float)(mxGetScalar(prhs[1]));
     maxIter=(int)round(mxGetScalar(prhs[2]));
    }
    
////////////////////////// First input.
    // First input should be x from (Ax=b), or the image.
    mxArray const * const image = prhs[0];
    mwSize const numDims = mxGetNumberOfDimensions(image);
    
    // Image should be dim 3
    if (numDims!=3){
        mexErrMsgIdAndTxt("err", "Image is not 3D");
    }
    // Now that input is ok, parse it to C data types.
    float const * const imgaux = static_cast<float const *>(mxGetData(image));
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    
    float *  img = (float*)malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float));
    // We need a float image, and, unfortunatedly, the only way of casting it is by value
    // Also, MATLAB is column mayor and C is row mayor! we need to deal with that
    for(int i=0;i<size_img[0];i++)
        for(int j=0;j<size_img[1];j++)
            for(int k=0;k<size_img[2];k++)
        img[i*size_img[0]*size_img[1]+j*size_img[1]+k]=(float)imgaux[k*size_img[1]*size_img[2]+j*size_img[2]+i];
    
    
    
    
      
    
    // Allocte output image
    float *  imgout = (float*)malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float));
    // call C function with the CUDA denoising
  
    const long imageSize[3]={size_img[0] ,size_img[1],size_img[2] };
    pocs_tv(img,imgout, alpha, imageSize, maxIter); 
    
    //prepareotputs
    plhs[0] = mxCreateNumericArray(3,size_img, mxSINGLE_CLASS, mxREAL);
    float *mxImgout =(float*) mxGetPr(plhs[0]);
    
    memcpy(mxImgout,imgout,size_img[0] *size_img[1] *size_img[2]*sizeof(float));
    //free memory
    free(img);
    free(imgout);
    
}