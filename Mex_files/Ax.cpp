/*
 * @file Ax.cpp
 * @brief MEX gateway for Cone Beam  projection.
 * author: Ander Biguri ander.biguri@gmail.com
 *
 */

#include "tmwtypes.h"
#include "mex.h"
#include "matrix.h"
#include "ray_interpolated_projection.hpp"
#include "Siddon_projection.hpp"
#include <string.h>
#include <time.h>
/**
 * MEX gateway
 */



void mexFunction(int  nlhs , mxArray *plhs[],
        int nrhs, mxArray const *prhs[])
{
//     clock_t begin, end;
//     begin = clock();
    
    
    
    char const * const errId = "CBCT:MEX:Ax:InvalidInput";
    char const * const errMsgInputs = "Invalid number of inputs to MEX file.";
    char const * const errMsgImg = "Invalid dimension size of image (x) to MEX file.";
    
    
    
    //Check amount of inputs
    if (nrhs<3 ||nrhs>4) {
        mexErrMsgIdAndTxt(errId, errMsgInputs);
    }
    //////////////////////////// 4rd argument is matched or un matched
    bool krylov_proj=false;
    if (nrhs==4){
        if ( mxIsChar(prhs[3]) != 1)
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input","4rd input shoudl be a string");
        
        /* copy the string data from prhs[0] into a C string input_ buf.    */
        char *krylov = mxArrayToString(prhs[3]);
        if (strcmp(krylov,"Krylov"))
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input","4rd input shoudl be Krylov");
        else
            krylov_proj=true;
    }
    ///////////////////////// 3rd argument: angle of projection.
    
    size_t mrows = mxGetM(prhs[2]);
    size_t nalpha = mxGetN(prhs[2]);
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
            !(mrows==1) ) {
        mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input",
                "Input alpha must be a noncomplex array.");
    }
    mxArray const * const ptralphas=prhs[2];
    
    
    double const * const alphas = static_cast<double const *>(mxGetData(ptralphas));
    ////////////////////////// First input.
    // First input should be x from (Ax=b), or the image.
    mxArray const * const image = prhs[0];
    mwSize const numDims = mxGetNumberOfDimensions(image);
    
    // Image should be dim 3
    if (numDims!=3){
        mexErrMsgIdAndTxt(errId, errMsgImg);
    }
    // Now that input is ok, parse it to C data types.
    double const * const imgaux = static_cast<double const *>(mxGetData(image));
    // We need a float image, and, unfortunatedly, the only way of casting it is by value
    const mwSize *size_img= mxGetDimensions(image); //get size of image
    
    float *  img = (float*)malloc(size_img[0] *size_img[1] *size_img[2]* sizeof(float));
    for (int i=0;i<size_img[0] *size_img[1] *size_img[2];i++)
        img[i]=(float)imgaux[i];
    
    
    ///////////////////// Second input argument,
    // Geometry structure that has all the needed geometric data.
    
    // IMPORTANT-> Make sure Matlab creates the struct in this order.
    const char *fieldnames[11];
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
    
    
    if(!mxIsStruct(prhs[1]))
        mexErrMsgIdAndTxt( "CBCT:MEX:Ax:InvalidInput",
                "Second input must be a structure.");
    
    int nfields = mxGetNumberOfFields(prhs[1]);
    if (nfields < 10 || nfields >11 )
        mexErrMsgIdAndTxt("CBCT:MEX:Ax:InvalidInput","there are missing or extra fields in the geometry");
    
    // Check that all names are good
    mxArray    *tmp;
    size_t ncols;
    bool offsetAllOrig=false;
     bool offsetAllDetec=false;
    for(int ifield=0; ifield<nfields; ifield++) {
        tmp=mxGetField(prhs[1],0,fieldnames[ifield]);
        if(tmp==NULL){
            mexPrintf("%s number: %d %s \n", "FIELD",ifield+1, fieldnames[ifield]);
            mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputname",
                    "Above field is missing. Check spelling. ");
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
            case 6:case 7:case 10:
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp);
                if (mrows!=1 || ncols!=1){
                    mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                    mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                            "Above field has wrong size! Should be 1x1!");
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
    double *  acc;
    int c;
    Geometry geo;
    geo.unitX=1;geo.unitY=1;geo.unitZ=1;
//     mexPrintf("%d \n",nfields);
    for(int ifield=0; ifield<nfields; ifield++) {
        tmp=mxGetField(prhs[1],0,fieldnames[ifield]);
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
                geo.sVoxelX=sVoxel[0];
                geo.sVoxelY=sVoxel[1];
                geo.sVoxelZ=sVoxel[2];
                break;
            case 2:
                dVoxel=(double *)mxGetData(tmp);
                geo.dVoxelX=dVoxel[0];
                geo.dVoxelY=dVoxel[1];
                geo.dVoxelZ=dVoxel[2];
                break;
            case 3:
                nDetec=(double *)mxGetData(tmp);
                geo.nDetecU=(int)nDetec[0];
                geo.nDetecV=(int)nDetec[1];
                break;
            case 4:
                sDetec=(double *)mxGetData(tmp);
                geo.sDetecU=sDetec[0];
                geo.sDetecV=sDetec[1];
                break;
            case 5:
                dDetec=(double *)mxGetData(tmp);
                geo.dDetecU=dDetec[0];
                geo.dDetecV=dDetec[1];
                break;
            case 6:
                DSD=(double *)mxGetData(tmp);
                geo.DSD=DSD[0];
                break;
            case 7:
                DSO=(double *)mxGetData(tmp);
                geo.DSO=DSO[0];
                
                break;
            case 8:
               
                geo.offOrigX=(double*)malloc(nalpha * sizeof(double));
                geo.offOrigY=(double*)malloc(nalpha * sizeof(double));
                geo.offOrigZ=(double*)malloc(nalpha * sizeof(double));
                
                offOrig=(double *)mxGetData(tmp);
                
                for (int i=0;i<nalpha;i++){
                    if (offsetAllOrig)
                        c=i;
                    else
                        c=0;
                    geo.offOrigX[i]=offOrig[0+3*c];
                    geo.offOrigY[i]=offOrig[1+3*c];
                    geo.offOrigZ[i]=offOrig[2+3*c];
                }
                break;
            case 9:
                geo.offDetecU=(double*)malloc(nalpha * sizeof(double));
                geo.offDetecV=(double*)malloc(nalpha * sizeof(double));
                
                offDetec=(double *)mxGetData(tmp);
                for (int i=0;i<nalpha;i++){
                    if (offsetAllDetec)
                        c=i;
                    else
                        c=0;
                    geo.offDetecU[i]=offDetec[0+2*c];
                    geo.offDetecV[i]=offDetec[1+2*c];
                }
                break;
            case 10:
                acc=(double*)mxGetData(tmp);
                geo.accuracy=acc[0];
                break;
            default:
                mexErrMsgIdAndTxt( "CBCT:MEX:Ax:unknown","This shoudl not happen. Weird");
                break;
                
        }
    }
    if (nfields==10)
        geo.accuracy=0.2;
    
    
    // Additional test
    if( (size_img[0]!=geo.nVoxelX)|(size_img[1]!=geo.nVoxelY)|(size_img[2]!=geo.nVoxelZ))
        mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input",
                "Image size and nVoxel are not same size.");
    
    
    size_t num_bytes = geo.nDetecU*geo.nDetecV * sizeof(double);
    
    double** result = (double**)malloc(nalpha * sizeof(double*));
    for (int i=0; i<nalpha ;i++)
        result[i]=(double*)malloc(geo.nDetecU*geo.nDetecV *sizeof(double));
    
    // call the real function
    
//     end = clock();
//     double time_input = (double)(end - begin) / CLOCKS_PER_SEC;
//     mexPrintf("Input time : %lf ms\n" ,time_input*1000);
    
//     begin = clock();
    if (krylov_proj){
        siddon_ray_projection(img,geo,result,alphas,nalpha);
    }
    else
    {
        projection(img,geo,result,alphas,nalpha);
    }
    
//     end = clock();
//     double time_code = (double)(end - begin) / CLOCKS_PER_SEC;
//     mexPrintf("Maths time total : %lf ms\n" ,time_code*1000);
    
    
    // Set outputs and exit
    
//     begin = clock();
    mwSize* outsize;
    outsize[0]=geo.nDetecV;
    outsize[1]=geo.nDetecU;
    outsize[2]= nalpha;
    
    plhs[0] = mxCreateNumericArray(3,outsize,mxDOUBLE_CLASS,mxREAL);
//     plhs[0] = mxCreateNumericMatrix(geo.nDetecU,geo.nDetecV, ncols, mxDOUBLE_CLASS, mxREAL);
    double *outProjections = mxGetPr(plhs[0]);
    
    
    
    for (int i=0; i<nalpha ;i++)
        memcpy(&outProjections[geo.nDetecU*geo.nDetecV*i], result[i], geo.nDetecU*geo.nDetecV*sizeof(double));
    
    
    
    for (int i=0; i<nalpha ;i++)
        free (result[i]);
    free(result);
    
    // Free image data
    free(img);
    
//     end = clock();
//     double time_out = (double)(end - begin) / CLOCKS_PER_SEC;
//     mexPrintf("Out time : %lf ms\n" ,time_out*1000);
    return;
    
    
    
}
