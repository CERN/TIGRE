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
void apply_filtration (const float* pfIn, size_t uiULen, size_t uiVLen, const float* pfFilter, float* pfOut, const GpuIds& gpuids) {
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
    const float fULInv = 1./uiULen;

    float* d_pfInOut = nullptr;
    cudaMalloc((void **)&d_pfInOut, uiLen * sizeof(float));
    cudaCheckErrors("apply_filtration fail cudaMalloc 1");
    cudaMemcpy(d_pfInOut, pfIn, uiLen* sizeof(float), cudaMemcpyHostToDevice); // Sync only. pfIn is not pinned.
    cudaCheckErrors("apply_filtration fail cudaMemcpy 1");

    size_t uiBufferSize = (uiULen+1)/2+1;    // Buffer size for R2C. See https://docs.nvidia.com/cuda/cufft/

	cufftHandle cudafftPlanFwd;
	cufftHandle cudafftPlanInv;
    const int iBatch = uiVLen;
	cufftResult_t fftresult;
    fftresult = cufftPlan1d(&cudafftPlanFwd, uiULen, CUFFT_R2C, iBatch);
    cudafftCheckError(fftresult, "apply_filtration fail cufftPlan1d 1");
    fftresult = cufftPlan1d(&cudafftPlanInv, uiULen, CUFFT_C2R, iBatch);
    cudafftCheckError(fftresult, "apply_filtration fail cufftPlan1d 2");

    float* d_pfFilter = nullptr;
    cudaMalloc((void **)&d_pfFilter, uiULen * sizeof(float));
    cudaCheckErrors("apply_filtration fail cudaMalloc 2");
    cudaMemcpy(d_pfFilter, pfFilter, uiULen * sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 2");

    cufftComplex* d_pcfWork = nullptr;
    cudaMalloc((void **)&d_pcfWork, uiBufferSize * uiVLen*sizeof(cufftComplex));
    cudaCheckErrors("apply_filtration fail cudaMalloc 3");

    {
        const int divU = 128;//PIXEL_SIZE_BLOCK;
        const int divV = 1;//PIXEL_SIZE_BLOCK;
        dim3 grid((uiULen+divU-1)/divU,(uiVLen+divV-1)/divV,1);
        dim3 block(divU,divV,1);
        cufftSetStream(cudafftPlanFwd, 0);
        cufftSetStream(cudafftPlanInv, 0);
        fftresult = cufftExecR2C (cudafftPlanFwd, d_pfInOut, d_pcfWork);
        cudafftCheckError(fftresult, "apply_filtration fail cufftExecR2C");
        ApplyFilter<<<grid, block>>>(d_pcfWork, uiBufferSize, uiVLen, d_pfFilter, fULInv);// Kernel d_pcfInOut = d_pcfInOut * pfFilter / uiULen
        fftresult = cufftExecC2R (cudafftPlanInv, d_pcfWork, d_pfInOut);
        cudafftCheckError(fftresult, "apply_filtration fail cufftExecC2R");
    }
    cudaMemcpy(pfOut, d_pfInOut, uiLen*sizeof(float), cudaMemcpyDeviceToHost);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 3");

    cudaFree(d_pcfWork); d_pcfWork = nullptr;
    cudaFree(d_pfInOut); d_pfInOut = nullptr;
    cudaFree(d_pfFilter); d_pfFilter = nullptr;
    cufftDestroy(cudafftPlanFwd);
    cufftDestroy(cudafftPlanInv);
}


//! Apply filter in the Fourier space
void apply_filtration2 (const float* pfInAll, size_t uiOffset, size_t uiULen, size_t uiBatch, const float* pfFilter, size_t uiFLen, float fScale, float* pfOut, const GpuIds& gpuids) {
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
    // USING THE FIRST GPU ONLY
    const float* pfIn = pfInAll+uiOffset;
    cudaSetDevice(gpuids[0]);
    cudaCheckErrors("apply_filtration fail cudaSetDevice");
    size_t uiPaddingLen = (uiFLen-uiULen) / 2;
    float* d_pfProjWide = nullptr;
    cudaMalloc((void**)&d_pfProjWide, uiFLen*uiBatch*sizeof(float));
    cudaCheckErrors("apply_filtration fail cudaMalloc wide");
    cudaMemset(d_pfProjWide, 0, uiFLen*uiBatch*sizeof(float));
    cudaCheckErrors("apply_filtration fail cudaMemset");
    cudaMemcpy2D(&d_pfProjWide[uiPaddingLen], uiFLen*sizeof(float), pfIn, uiULen*sizeof(float), uiULen*sizeof(float), uiBatch, cudaMemcpyHostToDevice);
    cudaCheckErrors("apply_filtration fail cudaMemcpy2D");

    const float fFLInv = 1./uiFLen;

    size_t uiBufferSize = (uiFLen+1)/2+1;    // Buffer size for R2C. See https://docs.nvidia.com/cuda/cufft/

	cufftHandle cudafftPlanFwd;
	cufftHandle cudafftPlanInv;
	cufftResult_t fftresult;
    fftresult = cufftPlan1d(&cudafftPlanFwd, uiFLen, CUFFT_R2C, uiBatch);
    cudafftCheckError(fftresult, "apply_filtration fail cufftPlan1d 1");
    fftresult = cufftPlan1d(&cudafftPlanInv, uiFLen, CUFFT_C2R, uiBatch);
    cudafftCheckError(fftresult, "apply_filtration fail cufftPlan1d 2");

    float* d_pfFilter = nullptr;
    cudaMalloc((void **)&d_pfFilter, uiFLen * sizeof(float));
    cudaCheckErrors("apply_filtration fail cudaMalloc 2");
    cudaMemcpy(d_pfFilter, pfFilter, uiFLen * sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 2");

    cufftComplex* d_pcfWork = nullptr;
    cudaMalloc((void **)&d_pcfWork, uiBufferSize * uiBatch*sizeof(cufftComplex));
    cudaCheckErrors("apply_filtration fail cudaMalloc 3");

    {
        const int divU = 128;//PIXEL_SIZE_BLOCK;
        const int divV = 1;//PIXEL_SIZE_BLOCK;
        dim3 grid((uiFLen+divU-1)/divU,(uiBatch+divV-1)/divV,1);
        dim3 block(divU,divV,1);
        cufftSetStream(cudafftPlanFwd, 0);
        cufftSetStream(cudafftPlanInv, 0);
        fftresult = cufftExecR2C (cudafftPlanFwd, d_pfProjWide, d_pcfWork);
        cudafftCheckError(fftresult, "apply_filtration fail cufftExecR2C");
        ApplyFilter<<<grid, block>>>(d_pcfWork, uiBufferSize, uiBatch, d_pfFilter, fFLInv*fScale);// Kernel d_pcfInOut = d_pcfInOut * pfFilter / uiFLen * 
        fftresult = cufftExecC2R (cudafftPlanInv, d_pcfWork, d_pfProjWide);
        cudafftCheckError(fftresult, "apply_filtration fail cufftExecC2R");
    }
    cudaMemcpy2D(pfOut, uiULen*sizeof(float), &d_pfProjWide[uiPaddingLen], uiFLen*sizeof(float), uiULen*sizeof(float), uiBatch, cudaMemcpyDeviceToHost);
    cudaCheckErrors("apply_filtration fail cudaMemcpy 3");

    cudaFree(d_pcfWork); d_pcfWork = nullptr;
    cudaFree(d_pfProjWide); d_pfProjWide = nullptr;
    cudaFree(d_pfFilter); d_pfFilter = nullptr;
    cufftDestroy(cudafftPlanFwd);
    cufftDestroy(cudafftPlanInv);
}