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
#include "cuda_to_hip.h"
#ifndef USE_HIP
#include <curand_kernel.h>
#include <curand.h>
#endif

#include "gpuUtils.hpp"
#include "RandomNumberGenerator.hpp"
#include "tigre_release_on_error.hpp"

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                TIGRE_MATLAB_ONLY(cudaDeviceReset());\
                tigreCudaError("RandomNumberGenerator:", msg, cudaGetErrorString(__err));\
        } \
} while (0)


__global__ void setup_kernel(curandState *state, unsigned long long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /* Each thread gets same seed, a different sequence number, no offset */
    curand_init(seed, idx, 0, &state[idx]);
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
// One random sequence per GPU thread, all from "seed": the same seed gives the same noise. TIGRE's callers draw a new
// seed for every call unless the user gives one (before, it was always 1234, so every call added the same noise).
static void noiseOnGpu(const float* pfIn, size_t uiLen, bool gaussian, float fGaussMu, float fGaussSigma, float* pfOut,
                       const GpuIds& gpuids, unsigned long long seed)
{
    if (gpuids.GetLength() > 0) cudaSetDevice(gpuids[0]);  // the first of the GPUs asked for (before: always the current one)
    float* d_pfIn = nullptr;
    float* d_pfOut = nullptr;
    curandState *curandStates = nullptr;
    // If a CUDA error raises (Python build), free what this call allocated; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    onError.add([&]{ cudaFree(d_pfIn); cudaFree(d_pfOut); cudaFree(curandStates); });
    const size_t bytes = uiLen*sizeof(float);
    tigre::checkMalloc(cudaMalloc((void **)&d_pfIn, bytes), "RandomNumberGenerator:", "Noise generation (addCTnoise)", bytes);
    tigre::checkMalloc(cudaMalloc((void **)&d_pfOut, bytes), "RandomNumberGenerator:", "Noise generation (addCTnoise)", bytes);
    cudaMemcpy(d_pfIn, pfIn, bytes, cudaMemcpyHostToDevice);
    cudaCheckErrors("fail cudaMemcpy 1");

    const int kiBlockDim = gaussian ? 64 : 1024;  // Threads per Block
    const int kiGridDim = 64;
    const size_t stateBytes = (size_t)kiGridDim * kiBlockDim * sizeof(curandState);
    tigre::checkMalloc(cudaMalloc((void **)&curandStates, stateBytes), "RandomNumberGenerator:", "Noise generation (addCTnoise)", stateBytes);
    setup_kernel<<<kiGridDim, kiBlockDim>>>(curandStates, seed);
    if (gaussian)
        GeneratePoissonAddGaussian<<<kiGridDim, kiBlockDim>>>(curandStates, d_pfIn, uiLen, fGaussMu, fGaussSigma, d_pfOut);
    else
        GeneratePoisson<<<kiGridDim, kiBlockDim>>>(curandStates, d_pfIn, uiLen, d_pfOut);
    cudaMemcpy(pfOut, d_pfOut, bytes, cudaMemcpyDeviceToHost);
    cudaCheckErrors("fail cudaMemcpy 2");

    onError.dismiss();
    cudaFree(d_pfIn);
    cudaFree(d_pfOut);
    cudaFree(curandStates);
}

void poisson_1d(const float* pfIn, size_t uiLen, float* pfOut, const GpuIds& gpuids, unsigned long long seed) {
    noiseOnGpu(pfIn, uiLen, false, 0.0f, 0.0f, pfOut, gpuids, seed);
}

void poisson_gaussian_1d(const float* pfIn,
                        size_t uiLen,
                        float fGaussMu,
                        float fGaussSigma,
                        float* pfOut,
                        GpuIds& gpuids,
                        unsigned long long seed)
{
    noiseOnGpu(pfIn, uiLen, true, fGaussMu, fGaussSigma, pfOut, gpuids, seed);
}
