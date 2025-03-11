/*-------------------------------------------------------------------------
 * CUDA function for optimized proton CT radiographies
 * The full method is described in Kaser et al.: Integration of proton imaging into the TIGRE toolbox (submitted to ZMP)
 * and based on the method of Collins-Fekete (https://doi.org/10.1088/0031-9155/61/23/8232)
 */
 
/*--------------------------------------------------------------------------
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

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <iostream>
#ifndef improvedForwardProjections_H
#define improvedForwardProjections_H
#define pi 3.14159265359
#define eps 1e-8
#define vecSizeCS 220
#define vecSizeOut 100
#define vecSizeIn 10
#define maxthreads 256
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>

void ParticleProjections(float* outProjection, float* posIn, float* posOut, float* dirIn, float* dirOut, float* p_wepl, \
        int numOfEntries, int detectSizeX, int detectSizeY, float* pixelSize, float detectDistIn, float detectDistOut, float ein, float* ch_param);

__device__ int calcIntercepts(float* InterceptsVec ,float*  a, float* b, \
                      float*  c, float* d, float* pos1, float pixelSize, bool* protFlag, int maxIntercep);

__device__ int SolvePolynomial(float*x, float a, float b, float c);

__device__ int MinMax(float* solutions, float a, float b, float c);

__device__ void SimpleSort(float* arr, int size_arr);

__global__ void ParticleKernel(float* dhist1, float* dhist2, float* devicePosIn, float* devicePosOut, float* devicedirIn, \
                               float* devicedirOut ,float* p_wepl,int* numOfEntries, int* detectSizeX, int *detectSizeY, \
                               float* pixelSize, float *detectDistIn, float *detectDistOut, float *ein, float *hull, float *reject);

__device__ int hullEntryExit(float* HullIntercept, float* position, float* direction, int in_or_out, float *hullparams, float detOff);

__device__ int calcInterceptsLinear(float* LinInterceptsVec, float* start, float* stop, float* direction, float pix, int maxIntercep, \
        bool* protFlag);

void ParticleProjectionsCone(float* outProjection, float* posIn, float* posOut, float* dirIn, float* dirOut, float* p_wepl, \
        int numOfEntries, int detectSizeX, int detectSizeY, float* pixelSize, float detectDistIn, float detectDistOut, float sourcePos, \
        float ein, float* ch_param);

__device__ int calcInterceptsCone(float* InterceptsVec ,float*  a, float* b, \
                      float*  c, float* d, float* pos1, float pixelSize, bool* protFlag, int maxIntercep, \
                      float sourcePos, float din, float dout);

__device__ int SolvePolynomialCone(float*x, float a, float b, float c);

__device__ void SimpleSortCone(float* arr, int size_arr);

__device__ int MinMaxCone(float* solutions, float a, float b, float c);

__global__ void ParticleKernelCone(float* dhist1, float* dhist2, float* devicePosIn, float* devicePosOut, float* devicedirIn, \
                               float* devicedirOut ,float* p_wepl,int* numOfEntries, int* detectSizeX, int *detectSizeY, \
                               float* pixelSize, float *detectDistIn, float *detectDistOut, float *ein, float *hull, float *reject, \
                               float* sourceDist);

__device__ int hullEntryExitCone(float* HullIntercept, float* position, float* direction, int in_or_out, float *hullparams, float detOff);

__device__ int calcInterceptsLinearCone(float* LinInterceptsVec, float* start, float* stop, float* direction, float pix, int maxIntercep, \
        bool* protFlag, float sourcePos);

#endif





















































































































































































