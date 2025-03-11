/*--------------------------------------------------------------------------
--------------------------------------------------------------------------
 This file is part of the TIGRE Toolbox
 
 Copyright (c) 2015, University of Bath and 
                     CERN-European Organization for Nuclear Research
                     All rights reserved.

 License:            Open Source under BSD. 
                     See the full license at
                     https://github.com/CERN/TIGRE/blob/master/LICENSE

 Contact:            tigre.toolbox@gmail.com
 Codes:              https://github.com/CERN/TIGRE/
 Coded by:           Stefanie Kaser, Benjamin Kirchmayer 
--------------------------------------------------------------------------*/

#include "mex.h"
#include "CUDA/improvedForwardProjections.hpp"
#include <stdexcept>
#include <iostream>
#include <cstring>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[]){

    if (nrhs =! 7){
        mexErrMsgIdAndTxt("CS Projections:", "Check Number of Input arguments!");
    }

    float *posIn, *posOut, *dirIn, *dirOut;
    float *Wepl, *pixelSize, *detectorDistanceIn, *detectorDistanceOut, *initEnergy;

    //Load parameters
    posIn = (float *)(mxGetPr(prhs[0]));
    posOut = (float *)mxGetPr(prhs[1]);
    dirIn = (float *)mxGetPr(prhs[2]);
    dirOut = (float *)mxGetPr(prhs[3]);
    Wepl = (float*) mxGetPr(prhs[4]);
    initEnergy = (float*) mxGetPr(prhs[5]);

    //Get Number of Protons contained in the root files
    int numOfProtons = (int) mxGetM(prhs[4]);

    mxArray * geometryMex=(mxArray*)prhs[6];

    const char *fieldnames_geo[7];
    fieldnames_geo[0] = "dDetector";
    fieldnames_geo[1] = "DSD";
    fieldnames_geo[2] = "DSID";
    fieldnames_geo[3] = "DSO";
    fieldnames_geo[4] = "hull";
    fieldnames_geo[5] = "sDetector";
    fieldnames_geo[6] = "mode";
    
    double * pix0, *dsd0, *dsid0, *hull0, *det0, *dso0;
    float pix[2], dsd, dsid, dso, hull[4], det[2];
    const char* mode;
    bool coneBeam = true;
    mxArray    *tmp;
    for (int ifield=0; ifield<7; ifield++){
        tmp=mxGetField(geometryMex,0,fieldnames_geo[ifield]);
        switch(ifield){
            case 0:
                pix0 =(double *)mxGetData(tmp);
                pix[0] = (float)pix0[0];
                pix[1] = (float)pix0[1];
                break;
            case 1:
                dsd0 =(double *)mxGetData(tmp);
                dsd = (float)dsd0[0];
                break;
            case 2:
                dsid0 =(double *)mxGetData(tmp);
                dsid = (float)dsid0[0];
                break;
            case 3:
                dso0 =(double *)mxGetData(tmp);
                dso = (float)dso0[0];
                break;
            case 4:
                hull0 =(double *)mxGetData(tmp);
                hull[0] = (float)hull0[0];
                hull[1] = (float)hull0[1];
                hull[2] = (float)hull0[2];
                hull[3] = (float)hull0[3];
                break;
            case 5:
                det0 =(double *)mxGetData(tmp);
                det[0] = (float)det0[0];
                det[1] = (float)det0[1];
                break;
            case 6:
                mode="";
                mode=mxArrayToString(tmp);
                if (!strcmp(mode,"parallel"))
                    coneBeam=false;
                break;
        } 
    }
    
    
    if (hull[3] == 0){std::cout << "Info: Calculation of optimized proton radiographies will be performed without object hull!" << std::endl;}
    
    if (hull[2] > 6.28318530717958648){std::cout << "Info: Hull rotation angle exceeds 2 Pi. Please check the input! Continuing with calculation..." << std::endl;}

    mwSize outSize[2];
    outSize[0] = int(det[1]/pix[1]);
    outSize[1] = int(det[0]/pix[0]);
    plhs[0] = mxCreateNumericArray(2, outSize, mxSINGLE_CLASS, mxREAL);
    float *outProjections = (float*)mxGetPr(plhs[0]);
   
    //For Calculation 2 historgrams are needed
    // 
    if(coneBeam == false){
        std::cout << "Info: Parallel geometry selected..." << std::endl;
        ParticleProjections(outProjections, posIn, posOut, dirIn, dirOut, Wepl, numOfProtons, int(det[0]/pix[0]), int(det[1]/pix[1]), pix, dsid-dso, dsd-dso, *initEnergy, hull);
    }
    else{
        std::cout << "Info: Cone beam geometry selected..." << std::endl;
        ParticleProjectionsCone(outProjections, posIn, posOut, dirIn, dirOut, Wepl, numOfProtons, int(det[0]/pix[0]), int(det[1]/pix[1]), pix, dsid-dso, dsd-dso, -1*dso, *initEnergy, hull);
    }

}
