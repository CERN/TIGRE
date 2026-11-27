/*-------------------------------------------------------------------------
 * gpuArray (Parallel Computing Toolbox) inputs and outputs for Ax_mex and Atb_mex.
 *
 * A gpuArray input runs on the GPU that holds it, straight from and to GPU
 * memory (the gpuData path of the CUDA projectors, as for PyTorch/CuPy arrays
 * in Python), and the result is a gpuArray on the same GPU.
 *
 * The toolbox's C API (mxGPUArray) is only compiled in with -DTIGRE_GPUARRAY,
 * which Compile.m adds when the toolbox is installed. Without it, a gpuArray
 * input gives an error that Ax.m and Atb.m catch to use a host copy instead.
 *--------------------------------------------------------------------------
 * This file is part of the TIGRE Toolbox
 *
 * Copyright (c) 2015, University of Bath and
 *                     CERN-European Organization for Nuclear Research
 *                     All rights reserved.
 *
 * License:            Open Source under BSD.
 *                     See the full license at
 *                     https://github.com/CERN/TIGRE/blob/master/LICENSE
 *
 * Contact:            tigre.toolbox@gmail.com
 * Codes:              https://github.com/CERN/TIGRE/
 *--------------------------------------------------------------------------*/
#ifndef TIGRE_GPUARRAY_HPP
#define TIGRE_GPUARRAY_HPP

#include <mex.h>
#include <CUDA/GpuIds.hpp>

inline bool tigreIsGpuArray(const mxArray* a) {
    return mxIsClass(a, "gpuArray");
}

#ifdef TIGRE_GPUARRAY

#include <gpu/mxGPUArray.h>
#include <cuda_runtime_api.h>

// A single, real gpuArray input, read only.
struct TigreGpuInput {
    mxGPUArray const* array = nullptr;
    float* data = nullptr;
    int device = -1;
    mwSize numDims = 0;
    mwSize dims[3] = {1, 1, 1};

    TigreGpuInput(const mxArray* a, const char* errId) {
        mxInitGPU();
        array = mxGPUCreateFromMxArray(a);
        if (mxGPUGetClassID(array) != mxSINGLE_CLASS || mxGPUGetComplexity(array) != mxREAL) {
            mxGPUDestroyGPUArray(array);
            array = nullptr;
            mexErrMsgIdAndTxt(errId, "gpuArray input should be real and single.");
        }
        data = static_cast<float*>(const_cast<void*>(mxGPUGetDataReadOnly(array)));
        numDims = mxGPUGetNumberOfDimensions(array);
        const mwSize* d = mxGPUGetDimensions(array);
        for (mwSize i = 0; i < numDims && i < 3; i++) dims[i] = d[i];
        mxFree(const_cast<mwSize*>(d));
        // Run on the GPU that holds the array (MATLAB's current GPU).
        cudaPointerAttributes attr;
        if (cudaPointerGetAttributes(&attr, data) != cudaSuccess || attr.type != cudaMemoryTypeDevice) {
            mxGPUDestroyGPUArray(array);
            array = nullptr;
            mexErrMsgIdAndTxt(errId, "Could not find the GPU that holds the gpuArray input.");
        }
        device = attr.device;
        // MATLAB queues gpuArray work asynchronously: make sure the input is written.
        cudaSetDevice(device);
        cudaDeviceSynchronize();
    }
    ~TigreGpuInput() {
        if (array) mxGPUDestroyGPUArray(array);
    }
    // Only the GPU that holds the array.
    void setGpuIds(GpuIds& gpuids) const {
        int dev = device;
        gpuids.SetIds(1, &dev);
    }
};

// A zero-filled single gpuArray output on the current GPU, returned as an mxArray.
struct TigreGpuOutput {
    mxGPUArray* array = nullptr;
    float* data = nullptr;

    TigreGpuOutput(mwSize numDims, const mwSize* dims) {
        array = mxGPUCreateGPUArray(numDims, dims, mxSINGLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES);
        data = static_cast<float*>(mxGPUGetData(array));
    }
    // Waits for TIGRE's kernels (they run on TIGRE's own streams) and hands the array to MATLAB.
    mxArray* release(int device) {
        cudaSetDevice(device);
        cudaDeviceSynchronize();
        mxArray* out = mxGPUCreateMxArrayOnGPU(array);
        mxGPUDestroyGPUArray(array);
        array = nullptr;
        return out;
    }
    ~TigreGpuOutput() {
        if (array) mxGPUDestroyGPUArray(array);
    }
};

#endif  // TIGRE_GPUARRAY
#endif  // TIGRE_GPUARRAY_HPP
