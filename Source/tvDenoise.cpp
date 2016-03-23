#include "tmwtypes.h"
#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "tvdenoising.hpp"
#include <string.h>
// #include <time.h>
/**
 * MEX gateway
 */
void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
    int maxIter;
    float lambda;
    if (nrhs==1){
        maxIter=100;
        lambda=15.0f;
    }
    if (nrhs==2){
       mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "Only 1 TV hyperparemter inputed");
    }
    if (nrhs>3){
       mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "Too many imput argumets");
    }
    if (nrhs==3){
     size_t mrows = mxGetM(prhs[1]);
     size_t ncols = mxGetN(prhs[1]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "TV parameters shoudl be 1x1");
     mrows = mxGetM(prhs[2]);
     ncols = mxGetN(prhs[2]);
     if (mrows!=1 || ncols !=1)
        mexErrMsgIdAndTxt("CBCT:CUDA:TVdenoising", "TV parameters shoudl be 1x1");
     lambda= (float)(mxGetScalar(prhs[1]));
     maxIter=(int)round(mxGetScalar(prhs[2]));
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
    double const * const imgaux = static_cast<double const *>(mxGetData(image));
    // We need a float image, and, unfortunatedly, the only way of casting it is by value
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    
    float *  img = (float*)malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float));
//     for (int i=0;i<size_img[0] *size_img[1] *size_img[2];i++)
//         img[i]=(float)imgaux[i];
    for(int i=0;i<size_img[0];i++)
        for(int j=0;j<size_img[1];j++)
            for(int k=0;k<size_img[2];k++)
               img[i*size_img[0]*size_img[1]+j*size_img[1]+k]=(float)imgaux[k*size_img[1]*size_img[2]+j*size_img[2]+i];
    //////////////
    
    // Allocte output image
    float *  imgout = (float*)malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float));
    // call C function with the CUDA denoising
    const float spacing[3]={1,1,1};
    const long imageSize[3]={size_img[0] ,size_img[1],size_img[2] };
   
    tvdenoising(img,imgout, lambda, spacing, imageSize, maxIter); 
    
    //prepareotputs
    plhs[0] = mxCreateNumericArray(3,size_img, mxSINGLE_CLASS, mxREAL);
    float *mxImgout =(float*) mxGetPr(plhs[0]);
    
    memcpy(mxImgout,imgout,size_img[0] *size_img[1] *size_img[2]*sizeof(float));
    //free memory
    free(img);
    free(imgout);
     

}