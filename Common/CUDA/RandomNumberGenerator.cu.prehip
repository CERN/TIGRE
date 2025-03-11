/*-------------------------------------------------------------------------
 *
 * CUDA functions for random number generator
 *
 * Adds noise of Poisson and normal distribution to the input.
 *
 * CODE by       Tomoyuki SADAKANE
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

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <curand.h>

#include "gpuUtils.hpp"
#include "RandomNumberGenerator.hpp"

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                cudaDeviceReset();\
                mexErrMsgIdAndTxt("RandomNumberGenerator:",cudaGetErrorString(__err));\
        } \
} while (0)


__global__ void setup_kernel(curandState *state) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /* Each thread gets same seed, a different sequence number, no offset */
    curand_init(1234, idx, 0, &state[idx]);
}

__global__ void GeneratePoisson(curandState *state, const float* pfIn, size_t uiLen, float* pfOut) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /* Copy state to local memory for efficiency */
    curandState localState = state[idx];
    int iIter = (uiLen + blockDim.x*gridDim.x - 1)/(blockDim.x*gridDim.x);
    for (int iI = 0; iI < iIter; ++iI) {
        size_t uiPos = (size_t)blockDim.x*gridDim.x*iI+idx;
        if (uiPos < uiLen) {
            /* Poisson */
            unsigned int uiPoisson = curand_poisson(&localState, pfIn[uiPos]);
            pfOut[uiPos] = (float)uiPoisson;
        }
    }
    /* Copy state back to global memory */
    state[idx] = localState;
}

__global__ void GeneratePoissonAddGaussian(curandState *state,
                        const float* pfIn,
                        size_t uiLen, 
                        float fGaussMu,
                        float fGaussSigma,
                        float* pfOut)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /* Copy state to local memory for efficiency */
    curandState localState = state[idx];
    int iIter = (uiLen + blockDim.x*gridDim.x - 1)/(blockDim.x*gridDim.x);
    for (int iI = 0; iI < iIter; ++iI) {
        size_t uiPos = (size_t)blockDim.x*gridDim.x*iI+idx;
        if (uiPos < uiLen) {
            /* Poisson */
            unsigned int uiPoisson = curand_poisson(&localState, pfIn[uiPos]);
            /* Gaussian */
            float fNormal = curand_normal(&localState) * fGaussSigma + fGaussMu;
            pfOut[uiPos] = fNormal + (float)uiPoisson;
        }
    }
    /* Copy state back to global memory */
    state[idx] = localState;
}


template<class T_value>
void GetMinMax(const T_value* pfIn, size_t uiLen, T_value& tvMin, T_value& tvMax) {
    tvMin = pfIn[0];
    tvMax = pfIn[0];
    T_value tvVal;
    for (int iI = 1; iI < uiLen; ++iI) {
        tvVal = pfIn[iI];
        if (tvMax < tvVal) { tvMax = tvVal; continue;}
        if (tvMin > tvVal) { tvMin = tvVal; continue;}
    }
}
void poisson_1d(const float* pfIn, size_t uiLen, float* pfOut, const GpuIds& gpuids) {
    // printf("poisson_1d(pfIn = %p, uiLen = %zd, pfOut = %p)\n", pfIn, uiLen, pfOut);
    float* d_pfIn = nullptr;
    float* d_pfOut = nullptr;
    cudaMalloc((void **)&d_pfIn, uiLen * sizeof(float));
    cudaCheckErrors("poisson_1d fail cudaMalloc 1");
    cudaMalloc((void **)&d_pfOut, uiLen * sizeof(float));
    cudaCheckErrors("poisson_1d fail cudaMalloc 2");
    cudaMemcpy(d_pfIn, pfIn, uiLen*sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckErrors("poisson_1d fail cudaMemcpy 1");

    // float fMin, fMax;
    // GetMinMax(pfIn, uiLen, fMin, fMax);
    // printf("fMin, fMax = %f, %f\n", fMin, fMax);
    curandState *curandStates = nullptr;
    const int kiBlockDim = 1024;  // Threads per Block
    const int kiGridDim = 64;//(uiLen+kiBlockDim-1)/kiBlockDim;
    cudaMalloc((void **)&curandStates, kiGridDim * kiBlockDim * sizeof(curandState));
    cudaCheckErrors("poisson_1d fail cudaMalloc 3");
    setup_kernel<<<kiGridDim, kiBlockDim>>>(curandStates);
    GeneratePoisson<<<kiGridDim, kiBlockDim>>>(curandStates, d_pfIn, uiLen, d_pfOut);
    cudaMemcpy(pfOut, d_pfOut, uiLen*sizeof(float), cudaMemcpyDeviceToHost);
    cudaCheckErrors("poisson_1d fail cudaMemcpy 2");
    // GetMinMax(pfOut, uiLen, fMin, fMax);
    // printf("fMin, fMax = %f, %f\n", fMin, fMax);
    
    cudaFree(d_pfIn); d_pfIn = nullptr;
    cudaFree(d_pfOut); d_pfOut = nullptr;
    cudaFree(curandStates); curandStates = nullptr;
}

void poisson_gaussian_1d(const float* pfIn,
                        size_t uiLen,
                        float fGaussMu,
                        float fGaussSigma,
                        float* pfOut,
                        GpuIds& gpuids)
{
    // printf("poisson_gaussian_1d(pfIn = %p, uiLen = %zd, fGaussMu = %+f, fGaussSigma = %f, pfOut = %p)\n", pfIn, uiLen, fGaussMu, fGaussSigma, pfOut);
    float* d_pfIn = nullptr;
    float* d_pfOut = nullptr;
    cudaMalloc((void **)&d_pfIn, uiLen * sizeof(float));
    cudaCheckErrors("poisson_gaussian_1d fail cudaMalloc 1");
    cudaMalloc((void **)&d_pfOut, uiLen * sizeof(float));
    cudaCheckErrors("poisson_gaussian_1d fail cudaMalloc 2");
    cudaMemcpy(d_pfIn, pfIn, uiLen*sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckErrors("poisson_gaussian_1d fail cudaMemcpy 1");

    // float fMin, fMax;
    // GetMinMax(pfIn, uiLen, fMin, fMax);
    // printf("fMin, fMax = %f, %f\n", fMin, fMax);
    curandState *curandStates = nullptr;
    const int kiBlockDim = 64;  // Threads per Block
    const int kiGridDim = 64;//(uiLen+kiBlockDim-1)/kiBlockDim;
    cudaMalloc((void **)&curandStates, kiGridDim * kiBlockDim * sizeof(curandState));
    cudaCheckErrors("poisson_gaussian_1d fail cudaMalloc 3");
    setup_kernel<<<kiGridDim, kiBlockDim>>>(curandStates);
    GeneratePoissonAddGaussian<<<kiGridDim, kiBlockDim>>>(curandStates, d_pfIn, uiLen, fGaussMu, fGaussSigma, d_pfOut);
    cudaMemcpy(pfOut, d_pfOut, uiLen*sizeof(float), cudaMemcpyDeviceToHost);
    cudaCheckErrors("poisson_gaussian_1d fail cudaMemcpy 2");
    // GetMinMax(pfOut, uiLen, fMin, fMax);
    // printf("fMin, fMax = %f, %f\n", fMin, fMax);


    cudaFree(d_pfIn); d_pfIn = nullptr;
    cudaFree(d_pfOut); d_pfOut = nullptr;
    cudaFree(curandStates); curandStates = nullptr;
}
