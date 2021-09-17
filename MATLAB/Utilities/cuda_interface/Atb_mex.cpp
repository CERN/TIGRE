
/*-------------------------------------------------------------------------
 *
 * MATLAB MEX gateway for backprojection
 *
 * This file gets the data from MATLAB, checks it for errors and then
 * parses it to C and calls the relevant C/CUDA fucntions.
 *
 * CODE by       Ander Biguri
 *
 * ---------------------------------------------------------------------------
 * ---------------------------------------------------------------------------
 * Copyright (c) 2015, University of Bath and CERN- European Organization for
 * Nuclear Research
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * ---------------------------------------------------------------------------
 *
 * Contact: tigre.toolbox@gmail.com
 * Codes  : https://github.com/CERN/TIGRE
 * ---------------------------------------------------------------------------
 */



#include <math.h>
#include <string.h>
#include <tmwtypes.h>
#include <mex.h>
#include <matrix.h>
#include <CUDA/voxel_backprojection.hpp>
#include <CUDA/voxel_backprojection2.hpp>
#include <CUDA/voxel_backprojection_parallel.hpp>
#include <CUDA/GpuIds.hpp>





/**
 * MEX gateway
 *
 * This function takes data from MATLAB and passes it to the MEX code.
 * It checks and casts the inputs and prepares teh outputs for MATLAB.
 *
 *
 */

void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[]){
    
    //Check amount of inputs
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("CBCT:MEX:Atb:InvalidInput", "Wrong number of inputs provided");
    }
    ////////////////////////////
    // 5th argument is array of GPU-IDs.
    GpuIds gpuids;
    {
        size_t iM = mxGetM(prhs[4]);
        if (iM != 1) {
            mexErrMsgIdAndTxt( "CBCT:MEX:Atb:unknown","5th parameter must be a row vector.");
            return;
        }
        size_t uiGpuCount = mxGetN(prhs[4]);
        if (uiGpuCount == 0) {
            mexErrMsgIdAndTxt( "CBCT:MEX:Atb:unknown","5th parameter must be a row vector.");
            return;
        }
        int* piGpuIds = (int*)mxGetData(prhs[4]);
        gpuids.SetIds(uiGpuCount, piGpuIds);
    }
    
    /*
     ** 4th argument is matched or un matched.
     */
    bool pseudo_matched=false; // Caled krylov, because I designed it for krylov case....
    /* copy the string data from prhs[0] into a C string input_ buf.    */
    char *krylov = mxArrayToString(prhs[3]);
    if (!strcmp(krylov,"matched")) // if its 0, they are the same
        pseudo_matched=true;

    /*
     ** Third argument: angle of projection.
     */
    size_t mrows,nangles;
    
    mrows = mxGetM(prhs[2]);
    nangles = mxGetN(prhs[2]);
    
    
    mxArray const * const ptrangles=prhs[2];
    
    
    double const * const anglesM= static_cast<double const *>(mxGetData(ptrangles));
    // just copy paste the data to a float array
    float  *  angles= (float*)malloc(nangles*mrows*sizeof(float));
    for (int i=0;i<nangles*mrows;i++){
        angles[i]=(float)anglesM[i];
    }
    /**
     *
     * First input: The projections
     */
    
    // First input should be b from (Ax=b) i.e. the projections.
    mxArray const * const image = prhs[0];                 // Get pointer of the data
    mwSize const numDims = mxGetNumberOfDimensions(image); // Get numer of Dimensions of input matrix.
    // Image should be dim 3
    if (!(numDims==3 && nangles>1) && !(numDims==2 && nangles==1) ){
        mexErrMsgIdAndTxt("CBCT:MEX:Atb:InvalidInput",  "Projection data is not the right size");
    }
    if( !mxIsSingle(prhs[0])) {
        mexErrMsgIdAndTxt("CBCT:MEX:Ax:InvalidInput",
                "Input image must be a single noncomplex array.");
    }
    // Now that input is ok, parse it to C data types.
    // NOTE: while Number of dimensions is the size of the matrix in Matlab, the data is 1D row-wise mayor.
    
    // We need a float image, and, unfortunatedly, the only way of casting it is by value
//     const mwSize *size_proj= mxGetDimensions(image); //get size of image
//     mrows = mxGetM(image);
//     nangles = mxGetN(image);
//     size_t size_proj2;
//     if (nangles==1)
//         size_proj2=1;
//     else
//         size_proj2=size_proj[2];
    
    
    float  *  projections= static_cast<float *>(mxGetData(image));
    
    
    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Second input: Geometry structure
     */
    mxArray * geometryMex=(mxArray*)prhs[1];
    
    // IMPORTANT-> Make sure Matlab creates the struct in this order.
    const char *fieldnames[14];
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
    fieldnames[13]= "rotDetector";
    // Make sure input is structure
    
    mxArray    *tmp;
    
    // Now we know that all the input struct is good! Parse it from mxArrays to
    // C structures that MEX can understand.
    
    double * nVoxel, *nDetec; //we need to cast these to int
    double * sVoxel, *dVoxel,*sDetec,*dDetec, *DSO, *DSD,*offOrig,*offDetec;
    double *acc, *COR,*rotDetector;
    const char* mode;
    bool coneBeam=true;
    Geometry geo;
    int c;
    geo.unitX=1;geo.unitY=1;geo.unitZ=1;
    for(int ifield=0; ifield<14; ifield++) {
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
                geo.DSD=(float*)malloc(nangles * sizeof(float));
                DSD=(double *)mxGetData(tmp);
                for (int i=0;i<nangles;i++){
                    geo.DSD[i]=(float)DSD[i];
                }
                break;
            case 7:
                geo.DSO=(float*)malloc(nangles * sizeof(float));
                DSO=(double *)mxGetData(tmp);
                for (int i=0;i<nangles;i++){
                    geo.DSO[i]=(float)DSO[i];
                }
                break;
            case 8:
                
                geo.offOrigX=(float*)malloc(nangles * sizeof(float));
                geo.offOrigY=(float*)malloc(nangles * sizeof(float));
                geo.offOrigZ=(float*)malloc(nangles * sizeof(float));
                
                offOrig=(double *)mxGetData(tmp);
                
                for (int i=0;i<nangles;i++){
                    c=i;
                    
                    geo.offOrigX[i]=(float)offOrig[0+3*c];
                    geo.offOrigY[i]=(float)offOrig[1+3*c];
                    geo.offOrigZ[i]=(float)offOrig[2+3*c];
                }
                break;
            case 9:
                geo.offDetecU=(float*)malloc(nangles * sizeof(float));
                geo.offDetecV=(float*)malloc(nangles * sizeof(float));
                
                offDetec=(double *)mxGetData(tmp);
                for (int i=0;i<nangles;i++){
                    c=i;
                    
                    geo.offDetecU[i]=(float)offDetec[0+2*c];
                    geo.offDetecV[i]=(float)offDetec[1+2*c];
                }
                break;
            case 10:
                acc=(double*)mxGetData(tmp);
                geo.accuracy=(float)acc[0];
                if (acc[0]<0.001)
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:Accuracy","Accuracy should be bigger than 0.001");
                
                break;
            case 11:
                mode="";
                mode=mxArrayToString(tmp);
                if (!strcmp(mode,"parallel"))
                    coneBeam=false;
                break;
            case 12:
                COR=(double*)mxGetData(tmp);
                geo.COR=(float*)malloc(nangles * sizeof(float));
                for (int i=0;i<nangles;i++){
                    
                    geo.COR[i]  = (float)COR[i];
                }
                break;
            case 13:
                geo.dRoll= (float*)malloc(nangles * sizeof(float));
                geo.dPitch=(float*)malloc(nangles * sizeof(float));
                geo.dYaw=  (float*)malloc(nangles * sizeof(float));
                
                rotDetector=(double *)mxGetData(tmp);
                
                for (int i=0;i<nangles;i++){
                    c=i;
                    geo.dYaw[i]  = (float)rotDetector[0+3*c];
                    geo.dPitch[i]= (float)rotDetector[1+3*c];
                    geo.dRoll[i] = (float)rotDetector[2+3*c];
                    
                }
                break;
            default:
                mexErrMsgIdAndTxt( "CBCT:MEX:Atb:unknown","This should not happen. Weird");
                break;
                
        }
    }
    
    /*
     * allocate memory for the output (No longer needed)
     */
    
//     float* result = (float*)malloc(geo.nVoxelX *geo.nVoxelY*geo.nVoxelZ*sizeof(float));
    
    
    /*
     * Call the CUDA kernel
     */
    mwSize imgsize[3];
    imgsize[0]=geo.nVoxelX;
    imgsize[1]=geo.nVoxelY;
    imgsize[2]=geo.nVoxelZ;
    plhs[0] = mxCreateNumericArray(3,imgsize, mxSINGLE_CLASS, mxREAL);
    float *result = (float *)mxGetPr(plhs[0]);
    
    
    
    // To know which backprojection to call, we also need to know if the rotation is the orthodox/standard circular
    // rotation around the Z axis, or of its something else. This is because the current backprojection for
    // curcular scans is optimized with a trick that assumes that the voxels in Z direction
    // on the image are aligned with the axis of rotation, to incearse memory latency.
    // This however does not apply in arbitrary axis of rotation cases.
    // TODO: test if we really need 2 different codes, or if running the accelerated code
    // with the worng assumptions will just result in a speed like the non-accelerated code,
    // without sacrificing speedup in the standard case.
    
   
    // Run the CUDA code.
    if (coneBeam){
        if (pseudo_matched){
            voxel_backprojection2(projections,geo,result,angles,nangles, gpuids);
        }else{
            voxel_backprojection(projections,geo,result,angles,nangles, gpuids);
        }
    }else{
        voxel_backprojection_parallel(projections,geo,result,angles,nangles, gpuids);
    }


    
    return;
}
