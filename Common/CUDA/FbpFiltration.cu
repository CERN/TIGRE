/*-------------------------------------------------------------------------
 *
 * CUDA functions for convolution
 *
 * Applies the convolution filter in the Fourier space.
 * The filter should be given in the Fourier transformed form.
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

#include "TIGRE_common.hpp"
#include "FbpFiltration.hpp"
#include <string>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("apply_filtration",cudaGetErrorString(__err));\
        } \
} while (0)

void cudafftCheckError(cufftResult_t fftResult, const std::string& rstrMsg) {
    std::string strError = "Unknown error";
    if (fftResult == CUFFT_SUCCESS){ return; }
    else if (fftResult == CUFFT_INVALID_PLAN  ) { strError = "The plan parameter is not a valid handle. Handle is not valid when the plan is locked."; }
    else if (fftResult == CUFFT_ALLOC_FAILED  ) { strError = "The allocation of GPU resources for the plan failed."; }
    else if (fftResult == CUFFT_INVALID_VALUE ) { strError = "One or more invalid parameters were passed to the API."; }
    else if (fftResult == CUFFT_INTERNAL_ERROR) { strError = "An internal driver error was detected."; }
    else if (fftResult == CUFFT_SETUP_FAILED  ) { strError = "The cuFFT library failed to initialize."; }
    else if (fftResult == CUFFT_INVALID_SIZE  ) { strError = "The nx or batch parameter is not a supported size."; }
    mexPrintf("%s \n", rstrMsg.c_str());
    mexErrMsgIdAndTxt("ApplyFiltration", strError.c_str());
}

__global__ void ApplyFilter(cufftComplex* pcfInOut, size_t uiULen, size_t uiVLen, float* pfFilter, float fULInv) {

    size_t uiU = threadIdx.x + blockIdx.x * blockDim.x;
    size_t uiV = threadIdx.y + blockIdx.y * blockDim.y;
    if (uiV >= uiVLen || uiU >= uiULen) {
        return;
    }
    pcfInOut[uiU+uiULen*uiV].x *= pfFilter[uiU]*fULInv;
    pcfInOut[uiU+uiULen*uiV].y *= pfFilter[uiU]*fULInv;
}

//! Apply filter in the Fourier space
void apply_filtration(const float* pfIn, size_t uiULen, size_t uiVLen, const float* pfFilter, float* pfOut, const GpuIds& gpuids){
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    cudaCheckErrors("Device query fail");
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("apply_filtration","There are no available device(s) that support CUDA\n");
    }
    //
    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning thrown)
    // Check the available devices, and if they are the same
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("apply_filtration","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    // USE THE FIRST GPU ONLY!!!!!!!!!!!!!!!!! 
    cudaSetDevice(gpuids[0]);

    const size_t uiLen = uiULen * uiVLen;
    cufftComplex* h_pcfInOut = (cufftComplex*)malloc(uiLen*sizeof(cufftComplex));

    if (!h_pcfInOut) {
        mexErrMsgIdAndTxt("ApplyFiltration", "apply_filtration fail cudaMallocHost 1");
    }
    for (int iV = 0; iV < uiVLen; ++iV) {
        for (int iU = 0; iU < uiULen; ++iU) {
            h_pcfInOut[iU+uiULen*iV] = cufftComplex{pfIn[iU+uiULen*iV], 0};
        }
    }

    const float fULInv = 1./uiULen;

	cufftHandle cudafftPlan;
    const int iBatch = uiVLen;
	cufftResult_t fftresult;
    fftresult = cufftPlan1d(&cudafftPlan, uiULen, CUFFT_C2C, iBatch);
    cudafftCheckError(fftresult, "apply_filtration fail cufftPlan1d 1");
    fftresult = cufftPlan1d(&cudafftPlan, uiULen, CUFFT_C2C, iBatch);
    cudafftCheckError(fftresult, "apply_filtration fail cufftPlan1d 2");

    float*        d_pfFilter = nullptr;
    cufftComplex* d_pcfInOut = nullptr;
    cudaMalloc((void **)&d_pcfInOut, uiLen * sizeof(cufftComplex));
    cudaCheckErrors("apply_filtration fail cudaMalloc 1");
    cudaMalloc((void **)&d_pfFilter, uiULen * sizeof(float));
    cudaCheckErrors("apply_filtration fail cudaMalloc 2");
    cudaMemcpy(d_pcfInOut, h_pcfInOut, uiLen* sizeof(cufftComplex), cudaMemcpyHostToDevice);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 1");
    cudaMemcpy(d_pfFilter, pfFilter, uiULen * sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 2");

    {
        const int divU = 128;//PIXEL_SIZE_BLOCK;
        const int divV = 1;//PIXEL_SIZE_BLOCK;
        dim3 grid((uiULen+divU-1)/divU,(uiVLen+divV-1)/divV,1);
        dim3 block(divU,divV,1);
        cufftSetStream(cudafftPlan, 0);
        fftresult = cufftExecC2C (cudafftPlan, d_pcfInOut, d_pcfInOut, CUFFT_FORWARD);
        cudafftCheckError(fftresult, "apply_filtration fail cufftExecC2C CUFFT_FORWARD");
        ApplyFilter<<<grid, block>>>(d_pcfInOut, uiULen, uiVLen, d_pfFilter, fULInv);// Kernel d_pcfInOut = d_pcfInOut * pfFilter / uiULen
        fftresult = cufftExecC2C (cudafftPlan, d_pcfInOut, d_pcfInOut, CUFFT_INVERSE);
        cudafftCheckError(fftresult, "apply_filtration fail cufftExecC2C CUFFT_INVERSE");
    }

    cudaMemcpy(h_pcfInOut, d_pcfInOut, uiLen*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 3");

    cudaFree(d_pcfInOut); d_pcfInOut = nullptr;
    cudaFree(d_pfFilter); d_pfFilter = nullptr;

    cufftSetStream(cudafftPlan, 0);
    cufftDestroy(cudafftPlan);

    for (int iV = 0; iV < uiVLen; ++iV) {
        for (int iU = 0; iU < uiULen; ++iU) {
            pfOut[iU+uiULen*iV] = h_pcfInOut[iU+uiULen*iV].x;
        }
    }
    free(h_pcfInOut); h_pcfInOut = nullptr;
}
