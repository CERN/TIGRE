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
#include "CUDA/improvedForwardProjections_cone.hpp"
#include <stdexcept>
#include<iostream>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[]){

    if (nrhs =! 13){
        mexErrMsgIdAndTxt("CS Projections:", "Check Number of Input arguments!");
    }

    float *posIn, *posOut, *dirIn, *dirOut, *sourcePos;
    float *Wepl, *pixelSize, *detectorDistanceIn, *detectorDistanceOut, *initEnergy, *hullParams;
    int detectorSizeX, detectorSizeY;

    //Load parameters
    posIn = (float *)(mxGetPr(prhs[0]));
    posOut = (float *)mxGetPr(prhs[1]);
    dirIn = (float *)mxGetPr(prhs[2]);
    dirOut = (float *)mxGetPr(prhs[3]);
    Wepl = (float*) mxGetPr(prhs[4]);
    pixelSize = (float*) mxGetPr(prhs[5]);
    detectorSizeX = (int) mxGetScalar(prhs[6]) / (*pixelSize);
    detectorSizeY = (int) mxGetScalar(prhs[7]) / (*pixelSize);
    sourcePos = (float*) mxGetPr(prhs[8]);
    detectorDistanceIn = (float*) mxGetPr(prhs[9]);
    detectorDistanceOut = (float*) mxGetPr(prhs[10]);
    initEnergy = (float*) mxGetPr(prhs[11]);
    hullParams = (float*) mxGetPr(prhs[12]);
    
    int hull_len = (int) mxGetM(prhs[12]);
    
    if (hull_len != 4){
        throw std::invalid_argument(" Please check convex hull parameters. 4 entries are needed: [a, b, alpha, h] for (x*cos(alpha)-z*sin(alpha))²/a² + (x*sin(alpha)+z*cos(alpha))²/b² = 1 with maximum height h. If you wish to perform the calculation without convex hull, make sure the last entry, h, equals zero." );
    }
    
    /* if (hullParams[4] > abs(*detectorDistanceIn) or hullParams[4] > abs(*detectorDistanceOut)){
        throw std::invalid_argument(" Radius of convex hull cannot exceed distance between detectors and origin!" );
    }*/
    
    if (hullParams[3] == 0){std::cout << "Info: Calculation of optimized proton radiographies will be performed without object hull!" << std::endl;}
    
    if (hullParams[2] > 6.28318530717958648){std::cout << "Info: Hull rotation angle exceeds 2 Pi. Please check the input! Continuing with calculation..." << std::endl;}
    
    if (*detectorDistanceIn > 0){
        throw std::invalid_argument(" Entry detector distance has to be < 0!" );
    }
    
    if (*detectorDistanceOut < 0){
        throw std::invalid_argument(" Exit detector distance has to be > 0!" );
    }
    
    if (*sourcePos > *detectorDistanceIn){
        throw std::invalid_argument(" Source position has to be < POS_DET_IN!" );
    }

    //Get Number of Protons contained in the root files
    int numOfProtons = (int) mxGetM(prhs[4]);
    mwSize outSize[2];
    outSize[0] = detectorSizeX;
    outSize[1] = detectorSizeY;
    plhs[0] = mxCreateNumericArray(2, outSize, mxSINGLE_CLASS, mxREAL);
    float *outProjections = (float*)mxGetPr(plhs[0]);
   
    //For Calculation 2 historgrams are needed

    ParticleProjections(outProjections, posIn, posOut, dirIn, dirOut, Wepl ,numOfProtons, detectorSizeX, detectorSizeY, *pixelSize, *detectorDistanceIn, *detectorDistanceOut, *sourcePos, *initEnergy, hullParams);

}