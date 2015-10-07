
/**
 * @file trilinear.cpp
 * @brief MEX gateway for trilinear texture memory interpolation.
 * author: Ander Biguri ander.biguri@gmail.com
 *
 */

#include "tmwtypes.h"
#include "mex.h"
#include "matrix.h"
#include "projection.hpp"
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
    if (nrhs!=3) {
        mexErrMsgIdAndTxt(errId, errMsgInputs);
    }
    
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
    
    // IMPORTANT-> MAke sure Matlab creates the struct in this order.
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
    size_t mrows,ncols;
    for(int ifield=0; ifield<nfields; ifield++) { 
        tmp=mxGetField(prhs[1],0,fieldnames[ifield]);
        if(tmp==NULL){
            mexPrintf("%s number: %d %s \n", "FIELD",ifield+1, fieldnames[ifield]);
                      mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputname",
                        "Above field is missing. Check spelling. ");
        }
        switch(ifield){
            
            // cases where we want 3 input arrays.
            case 0:case 1:case 2:case 8:
                mrows = mxGetM(tmp);
                ncols = mxGetN(tmp); 
                if (mrows!=3 || ncols!=1){
                      mexPrintf("%s %s \n", "FIELD: ", fieldnames[ifield]);
                      mexErrMsgIdAndTxt( "CBCT:MEX:Ax:inputsize",
                        "Above field has wrong size! Should be 3x1!");
                }
                
                break;
            // this ones should be 2x1
            case 3:case 4:case 5:case 9:
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
    double * sVoxel, *dVoxel,*sDetec,*dDetec, *DSO, *DSD,*offOrig,*offDetec;
    double *  acc;
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
             case 8:
                 offOrig=(double *)mxGetData(tmp);
                 geo.offOrigX=offOrig[0];
                 geo.offOrigY=offOrig[1];
                 geo.offOrigZ=offOrig[2];
                 break;
             case 9:
                  offDetec=(double *)mxGetData(tmp);
                 geo.offDetecU=offDetec[0];
                 geo.offDetecV=offDetec[1];

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
    // 3rd argument: angle of projection.
   
    mrows = mxGetM(prhs[2]);
    ncols = mxGetN(prhs[2]);
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        !(mrows==1) ) {
      mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input",
        "Input alpha must be a noncomplex array.");
    }
    mxArray const * const ptralphas=prhs[2];
    
    double const * const alphas = static_cast<double const *>(mxGetData(ptralphas));

    // Additional test
    if( (size_img[0]!=geo.nVoxelX)|(size_img[1]!=geo.nVoxelY)|(size_img[2]!=geo.nVoxelZ))
     mexErrMsgIdAndTxt( "CBCT:MEX:Ax:input",
        "Image size and nVoxel are not same size.");
    

    size_t num_bytes = geo.nDetecU*geo.nDetecV * sizeof(double);
    
    double** result = (double**)malloc(ncols * sizeof(double*));
    for (int i=0; i<ncols ;i++)
        result[i]=(double*)malloc(geo.nDetecU*geo.nDetecV *sizeof(double));

    // call the real function
    
//     end = clock();
//     double time_input = (double)(end - begin) / CLOCKS_PER_SEC;
//     mexPrintf("Input time : %lf ms\n" ,time_input*1000);
    
//     begin = clock();

    projection(img,geo,result,alphas,ncols);
    
//     end = clock();
//     double time_code = (double)(end - begin) / CLOCKS_PER_SEC;
//     mexPrintf("Maths time total : %lf ms\n" ,time_code*1000);
    
    
    // Set outputs and exit
    
//     begin = clock();
    mwSize* outsize;
    outsize[0]=geo.nDetecV;
    outsize[1]=geo.nDetecU;
    outsize[2]= ncols;

    plhs[0] = mxCreateNumericArray(3,outsize,mxDOUBLE_CLASS,mxREAL);
//     plhs[0] = mxCreateNumericMatrix(geo.nDetecU,geo.nDetecV, ncols, mxDOUBLE_CLASS, mxREAL);
    double *outProjections = mxGetPr(plhs[0]);
    
 

    for (int i=0; i<ncols ;i++)
        memcpy(&outProjections[geo.nDetecU*geo.nDetecV*i], result[i], geo.nDetecU*geo.nDetecV*sizeof(double));
    

    
    for (int i=0; i<ncols ;i++)
        free (result[i]);
    free(result);
    
    // Free image data
    free(img);
    
//     end = clock();
//     double time_out = (double)(end - begin) / CLOCKS_PER_SEC;
//     mexPrintf("Out time : %lf ms\n" ,time_out*1000);
    return;
    

    
}
