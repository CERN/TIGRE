/*-------------------------------------------------------------------------
 *
 * MATLAB MEX gateway for projection
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
#include "matrix.h"
#include "ray_interpolated_projection.hpp"
#include "ray_interpolated_projection_parallel.hpp"
#include "Siddon_projection.hpp"
#include "Siddon_projection_parallel.hpp"
#include <string.h>

/**
 * MEX gateway
 */



void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
//     clock_t begin, end;
//     begin = clock();
    
    
    //Check amount of inputs
    if (nrhs<3 || nrhs>4) {
        mexErrMsgIdAndTxt("CBCT:MEX:Ax:InvalidInput", "Invalid number of inputs to MEX file.");
    }
    ////////////////////////////
    //4rd argument is interpolated or ray-voxel
    bool interpolated=false;
    if (nrhs==4){
        if ( mxIsChar(prhs[3]) != 1)
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:InvalidInput","4rd input shoudl be a string");
        
        /* copy the string data from prhs[0] into a C string input_ buf.    */
        char *krylov = mxArrayToString(prhs[3]);
        if (strcmp(krylov,"interpolated") && strcmp(krylov,"ray-voxel"))
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:InvalidInput","4rd input shoudl be either 'interpolated' or 'ray-voxel'");
        else
            // If its not ray-voxel, its "interpolated"
            if (~strcmp(krylov,"ray-voxel"))
                interpolated=true;
    }
    ///////////////////////// 3rd argument: angle of projection.
    
    size_t mrows = mxGetM(prhs[2]);
    size_t nalpha = mxGetN(prhs[2]);
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
            !(mrows==1) ) {
        mexErrMsgIdAndTxt("CBCT:MEX:Ax:InvalidInput",
                "Input alpha must be a noncomplex array.");
    }
    mxArray const * const ptralphas=prhs[2];
    
    
    double const * const alphasM= static_cast<double const *>(mxGetData(ptralphas));
    // just copy paste the data to a float array
    float  *  alphas= (float*)malloc(nalpha*sizeof(float));
    for (int i=0;i<nalpha;i++)
        alphas[i]=(float)alphasM[i];
    
    
    ////////////////////////// First input.
    // First input should be x from (Ax=b), or the image.
    mxArray const * const image = prhs[0];
    mwSize const numDims = mxGetNumberOfDimensions(image);
    
    // Image should be dim 3
    if (numDims!=3){
        mexErrMsgIdAndTxt( "CBCT:MEX:Ax:InvalidInput", "Invalid dimension size of image (x) to MEX file.");
    }
     if( !mxIsSingle(prhs[0])) {
       mexErrMsgIdAndTxt("CBCT:MEX:Ax:InvalidInput",
                "Input image must be a single noncomplex array.");
     }
    // Now that input is ok, parse it to C data types.
    float const * const imgaux = static_cast<float const *>(mxGetData(image));
    // We need a float image, and, unfortunatedly, the only way of casting it is by value
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    
    float *  img = (float*)malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float));
    for (unsigned long long i=0;i<size_img[0] *size_img[1] *size_img[2];i++)
        img[i]=(float)imgaux[i];
    
    ///////////////////// Second input argument,
    // Geometry structure that has all the needed geometric data.
    
    
    mxArray * geometryMex=(mxArray*)prhs[1];

    // IMPORTANT-> Make sure Matlab creates the struct in this order.
    const char *fieldnames[13];
    fieldnames[0] = "nVoxel";
    fieldnames[1] = "sVoxel";
    fieldnames[2] = "dVoxel";
    fieldnames[3] = "nDetector";
    fieldnames[4] = "sDetector";
    fieldnames[5] = "dDetector";
    fieldnames[6] = "DSD";
    fieldnames[7] = "DSO";
    fieldnames[8] = "offOrigin";
    fieldnames[9] = "offDetector";
    fieldnames[10]= "accuracy";
    fieldnames[11]= "mode";
    fieldnames[12]= "COR";

    
    if(!mxIsStruct(geometryMex))
        mexErrMsgIdAndTxt( "CBCT:MEX:Ax:InvalidInput",
                "Second input must be a structure.");
    int nfields = mxGetNumberOfFields(geometryMex);
    if (nfields < 10 || nfields >13 ){
         
        mexErrMsgIdAndTxt("CBCT:MEX:Ax:InvalidInput","there are missing or extra fields in the geometry");
    }
    // Check that all names are good
    mxArray    *tmp;
    size_t ncols;
    bool offsetAllOrig=false;
    bool offsetAllDetec=false;
    for(int ifield=0; ifield<13; ifield++) {
        tmp=mxGetField(geometryMex,0,fieldnames[ifield]);
        if(tmp==NULL){
           //tofix
            continue;
        }
        switch(ifield){
            
            // cases where we want 3 input arrays.
            case 0:case 1:case 2:
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp);
                if (mrows!=3 || ncols!=1){
                    
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexPrintf("%d x %d \n", "FIELD: ", (int)mrows,(int)ncols);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field has wrong size! Should be 3x1!");
                }
                
                break;
                //this one can be either 3x1 or 3xNangles
            case 8:
   
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp);

                if (mrows!=3 || ( ncols!=1&& ncols!=nalpha) ){
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexPrintf("%ld x %ld \n", "FIELD: ", (long int)mrows,(long int)ncols);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field has wrong size! Should be 3x1 or 3xlength(angles)!");
                    
                }
               
                if (ncols==nalpha)
                    offsetAllOrig=true;
                break;
                
            case 9:
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp);
                if (mrows!=2 || ( ncols!=1&& ncols!=nalpha)){
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field has wrong size! Should be 2x1 or 2xlength(angles)!");
                   
                }
                
                if (ncols==nalpha)
                    offsetAllDetec=true;
                break;
                // this ones should be 2x1
            case 3:case 4:case 5:
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp);
                if (mrows!=2 || ncols!=1){
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field has wrong size! Should be 2x1!");
                }
                break;
                // this ones should be 1x1
            case 6:case 7:case 10:case 12:
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp);
                if (mrows!=1 || ncols!=1){
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field has wrong size! Should be 1x1!");
                }
                
                break;
            case 11:
                if (!mxIsChar(tmp)){
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field is not string!");
                }
                break;
            default:
                mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input",
                        "something wrong happened. Ensure Geometric struct has correct amount of inputs.");
        }
        
    }
    // Now we know that all the input struct is good! Parse it from mxArrays to
    // C structures that MEX can understand.
    double * nVoxel, *nDetec; //we need to cast these to int
    double * sVoxel, *dVoxel,*sDetec,*dDetec, *DSO, *DSD;
    double *offOrig,*offDetec;
    double *  acc, *COR;
    const char* mode;
    int c;
    Geometry geo;
    geo.unitX=1;geo.unitY=1;geo.unitZ=1;
    bool coneBeam=true;
//     mexPrintf("%d \n",nfields);
    for(int ifield=0; ifield<13; ifield++) {
        tmp=mxGetField(geometryMex,0,fieldnames[ifield]);
         if(tmp==NULL){
           //tofix
            continue;
        }
        switch(ifield){
            case 0:
                nVoxel=(double *)mxGetData(tmp);
                // copy data to MEX memory
                geo.nVoxelX=(int)nVoxel[0];
                geo.nVoxelY=(int)nVoxel[1];
                geo.nVoxelZ=(int)nVoxel[2];
                break;
            case 1:
                sVoxel=(double *)mxGetData(tmp);
                geo.sVoxelX=(float)sVoxel[0];
                geo.sVoxelY=(float)sVoxel[1];
                geo.sVoxelZ=(float)sVoxel[2];
                break;
            case 2:
                dVoxel=(double *)mxGetData(tmp);
                geo.dVoxelX=(float)dVoxel[0];
                geo.dVoxelY=(float)dVoxel[1];
                geo.dVoxelZ=(float)dVoxel[2];
                break;
            case 3:
                nDetec=(double *)mxGetData(tmp);
                geo.nDetecU=(int)nDetec[0];
                geo.nDetecV=(int)nDetec[1];
                break;
            case 4:
                sDetec=(double *)mxGetData(tmp);
                geo.sDetecU=(float)sDetec[0];
                geo.sDetecV=(float)sDetec[1];
                break;
            case 5:
                dDetec=(double *)mxGetData(tmp);
                geo.dDetecU=(float)dDetec[0];
                geo.dDetecV=(float)dDetec[1];
                break;
            case 6:
                DSD=(double *)mxGetData(tmp);
                geo.DSD=(float)DSD[0];
                break;
            case 7:
                DSO=(double *)mxGetData(tmp);
                geo.DSO=(float)DSO[0];
                
                break;
            case 8:
               
                geo.offOrigX=(float*)malloc(nalpha * sizeof(float));
                geo.offOrigY=(float*)malloc(nalpha * sizeof(float));
                geo.offOrigZ=(float*)malloc(nalpha * sizeof(float));
                
                offOrig=(double *)mxGetData(tmp);
                
                for (int i=0;i<nalpha;i++){
                    if (offsetAllOrig)
                        c=i;
                    else
                        c=0;
                    geo.offOrigX[i]=(float)offOrig[0+3*c];
                    geo.offOrigY[i]=(float)offOrig[1+3*c];
                    geo.offOrigZ[i]=(float)offOrig[2+3*c];
                }
                break;
            case 9:
                geo.offDetecU=(float*)malloc(nalpha * sizeof(float));
                geo.offDetecV=(float*)malloc(nalpha * sizeof(float));
                
                offDetec=(double *)mxGetData(tmp);
                for (int i=0;i<nalpha;i++){
                    if (offsetAllDetec)
                        c=i;
                    else
                        c=0;
                    geo.offDetecU[i]=(float)offDetec[0+2*c];
                    geo.offDetecV[i]=(float)offDetec[1+2*c];
                }
                break;
            case 10:
                acc=(double*)mxGetData(tmp);
                if (acc[0]<0.001)
                   mexErrMsgIdAndTxt( "CBCT:MEX:Ax:Accuracy","Accuracy should be bigger than 0");
                   
                geo.accuracy=(float)acc[0];
                break;
            case 11:
                mode="";
                mode=mxArrayToString(tmp);
                if (!strcmp(mode,"parallel"))
                    coneBeam=false;
                else if (strcmp(mode,"cone"))
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:Mode","Unkown mode. Should be parallel or cone");
                break; 
             case 12:
                COR=(double*)mxGetData(tmp);
                geo.COR=(float)COR[0];
                break;
            default:
                mexErrMsgIdAndTxt( "CBCT:MEX:Ax:unknown","This shoudl not happen. Weird");
                break;
                
        }
    }
    tmp=mxGetField(geometryMex,0,fieldnames[10]);
    if (tmp==NULL)
        geo.accuracy=0.5;
    // Geometry
    tmp=mxGetField(geometryMex,0,fieldnames[11]);
    if (tmp==NULL)
        coneBeam=true;
    // COR
    tmp=mxGetField(geometryMex,0,fieldnames[12]);
    if (tmp==NULL)
        geo.COR=0.0;
    // Additional test
    if( (size_img[0]!=geo.nVoxelX)|(size_img[1]!=geo.nVoxelY)|(size_img[2]!=geo.nVoxelZ))
        mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input",
                "Image size and nVoxel are not same size.");
    
    size_t num_bytes = geo.nDetecU*geo.nDetecV * sizeof(float);
    
    float** result = (float**)malloc(nalpha * sizeof(float*));
    for (int i=0; i<nalpha ;i++)
        result[i]=(float*)malloc(geo.nDetecU*geo.nDetecV *sizeof(float));
    
    // call the real function
    
    if (coneBeam){
        if (interpolated){
            siddon_ray_projection(img,geo,result,alphas,nalpha);
        }else{
            interpolation_projection(img,geo,result,alphas,nalpha);
        }
    }else{
        if (interpolated){
//             mexErrMsgIdAndTxt( "CBCT:MEX:Ax:debug",
//                             "ray-voxel intersection is still unavailable for parallel beam, as there are some bugs on it.");
            siddon_ray_projection_parallel(img,geo,result,alphas,nalpha);
        }else{
            interpolation_projection_parallel(img,geo,result,alphas,nalpha);
        }
    }
    // Set outputs and exit
    
    mwSize outsize[3];
    outsize[0]=geo.nDetecV;
    outsize[1]=geo.nDetecU;
    outsize[2]= nalpha;
    
    plhs[0] = mxCreateNumericArray(3,outsize,mxSINGLE_CLASS,mxREAL);
    float *outProjections = (float*)mxGetPr(plhs[0]);
    for (int i=0; i<nalpha ;i++)
        for (unsigned long j=0; j<geo.nDetecU*geo.nDetecV;j++)
            outProjections[geo.nDetecU*geo.nDetecV*i+j]=(float)result[i][j];
            
    
    for (unsigned int i=0; i<nalpha ;i++)
        free (result[i]);
    free(result);
    
    // Free image data
    free(img);
    return; 
    
}
