/*-------------------------------------------------------------------------
 *
 * cuda_to_hip.h: CUDA -> HIP/ROCm compatibility shim for the Python/CUDA
 * backend of TIGRE.
 *
 * When USE_HIP is defined (the BUILD_WITH_HIP branch in setup.py) this header
 * maps the small set of CUDA runtime / cuRAND symbols the GPU sources use onto
 * their hip* equivalents and provides a texture fetch (tex3D_TIGRE) that uses
 * hardware trilinear filtering where available and a software lerp where it is
 * not. On NVIDIA the header is a thin passthrough to the CUDA runtime, so the
 * CUDA build is unchanged.
 *
 * fp32 element-read (cudaReadModeElementType) hardware linear filtering is not
 * available on every ROCm / hardware combination: some reject creation of a
 * cudaFilterModeLinear texture, others accept it but silently point-sample.
 * tigre_hw_linear_supported() (in gpuUtils.cu) decides at runtime with a
 * verified self-test -- it creates a tiny Linear ramp texture and confirms a
 * known sample actually interpolates -- caches the verdict for the process and
 * logs it once. On CUDA it returns true without testing.
 *
 * The backprojection kernels and the interpolated forward projector pick their
 * texture filterMode from that verdict (Linear when supported, else Point), and
 * tex3D_TIGRE mirrors it through the per-TU device flag tigre_hw_linear_filter:
 * when supported it forwards to the hardware tex3D<float>; otherwise it point-
 * samples the 8 neighbours through the point-filtered texture (border addressing
 * returns 0 outside the array, matching cudaAddressModeBorder) and lerps with
 * CUDA's unnormalized -0.5 texel-center convention. The texture's creation mode
 * and the sampling path are both driven from the single cached verdict so they
 * never disagree. tigre_sync_hw_linear() copies the verdict into this TU's
 * device flag and must be called once after the texture is created.
 *
 * CODE by       Jeff Daily <jeff.daily@amd.com>
 * ---------------------------------------------------------------------------
 * ---------------------------------------------------------------------------
 * Copyright (c) 2026, Advanced Micro Devices, Inc.
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
#ifndef CUDA_TO_HIP_H
#define CUDA_TO_HIP_H

#if defined(USE_HIP)

/* Host string/alloc functions must resolve to the C library, not the HIP
 * __device__ overloads, when a .cu is compiled as HIP host+device. */
#include <cstring>
#include <cstdlib>

#if defined(__HIP_DEVICE_COMPILE__) || defined(__HIPCC__)
#include <hip/hip_runtime.h>
#else
#include <hip/hip_runtime_api.h>
#endif

/* ---- runtime types / enums ---- */
#define cudaArray                       hipArray
#define cudaArray_t                     hipArray_t
#define cudaExtent                      hipExtent
#define cudaPitchedPtr                  hipPitchedPtr
#define cudaPos                         hipPos
#define cudaMemcpy3DParms               hipMemcpy3DParms
#define cudaChannelFormatDesc           hipChannelFormatDesc
#define cudaResourceDesc                hipResourceDesc
#define cudaTextureDesc                 hipTextureDesc
#define cudaTextureObject_t             hipTextureObject_t
#define cudaStream_t                    hipStream_t
#define cudaEvent_t                     hipEvent_t
#define cudaError_t                     hipError_t
#define cudaDeviceProp                  hipDeviceProp_t

#define cudaResourceTypeArray           hipResourceTypeArray
#define cudaReadModeElementType         hipReadModeElementType
#define cudaFilterModePoint             hipFilterModePoint
#define cudaFilterModeLinear            hipFilterModeLinear
#define cudaAddressModeBorder           hipAddressModeBorder
#define cudaAddressModeClamp            hipAddressModeClamp
#define cudaAddressModeWrap             hipAddressModeWrap

#define cudaMemcpyHostToDevice          hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost          hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice        hipMemcpyDeviceToDevice
#define cudaMemcpyHostToHost            hipMemcpyHostToHost
#define cudaMemcpyDefault               hipMemcpyDefault

#define cudaSuccess                     hipSuccess
#define cudaErrorMemoryAllocation       hipErrorOutOfMemory

#define cudaHostRegisterPortable        hipHostRegisterPortable
#define cudaDevAttrHostRegisterSupported hipDeviceAttributeHostRegisterSupported

/* ---- runtime API ---- */
#define cudaSetDevice                   hipSetDevice
#define cudaGetDevice                   hipGetDevice
#define cudaGetDeviceCount              hipGetDeviceCount
#define cudaGetDeviceProperties         hipGetDeviceProperties
#define cudaDeviceGetAttribute          hipDeviceGetAttribute
#define cudaDeviceSynchronize           hipDeviceSynchronize
#define cudaDeviceReset                 hipDeviceReset
#define cudaGetLastError                hipGetLastError
#define cudaPeekAtLastError             hipPeekAtLastError
#define cudaGetErrorString              hipGetErrorString
#define cudaMemGetInfo                  hipMemGetInfo

#define cudaMalloc                      hipMalloc
#define cudaFree                        hipFree
#define cudaMallocHost                  hipHostMalloc
#define cudaFreeHost                    hipHostFree
#define cudaHostRegister                hipHostRegister
#define cudaHostUnregister              hipHostUnregister
#define cudaMemset                      hipMemset
#define cudaMemsetAsync                 hipMemsetAsync
#define cudaMemcpy                      hipMemcpy
#define cudaMemcpyAsync                 hipMemcpyAsync
#define cudaMemcpyToSymbol              hipMemcpyToSymbol
#define cudaMemcpyToSymbolAsync         hipMemcpyToSymbolAsync

/* hipMalloc3DArray requires the flags argument that cudaMalloc3DArray
 * defaults to 0. */
#define cudaMalloc3DArray(arr, desc, ext) hipMalloc3DArray((arr), (desc), (ext), 0)
#define cudaFreeArray                   hipFreeArray
#define cudaMemcpy3D                    hipMemcpy3D
#define cudaMemcpy3DAsync               hipMemcpy3DAsync
#define cudaCreateChannelDesc           hipCreateChannelDesc
#define cudaCreateTextureObject         hipCreateTextureObject
#define cudaDestroyTextureObject        hipDestroyTextureObject

#define make_cudaExtent                 make_hipExtent
#define make_cudaPitchedPtr             make_hipPitchedPtr
#define make_cudaPos                    make_hipPos

#define cudaStreamCreate                hipStreamCreate
#define cudaStreamDestroy               hipStreamDestroy
#define cudaStreamSynchronize           hipStreamSynchronize
#define cudaStreamCreateWithFlags       hipStreamCreateWithFlags
#define cudaStreamNonBlocking           hipStreamNonBlocking

#define cudaEventCreate                 hipEventCreate
#define cudaEventDestroy                hipEventDestroy
#define cudaEventRecord                 hipEventRecord
#define cudaEventSynchronize            hipEventSynchronize
#define cudaEventElapsedTime            hipEventElapsedTime

/* HIP provides the round-to-nearest single-precision intrinsics but not the
 * round-toward-negative-infinity (_rd) variants. The reconstruction is
 * iterative and graded by nRMSE/adjointness tolerances, not bit-exact output,
 * so the 1-ULP rounding-mode difference is immaterial here. */
#define __fsqrt_rd                      __fsqrt_rn
#define __frcp_rd                       __frcp_rn

/* ---- cuRAND -> hipRAND ---- */
#if defined(__HIP_DEVICE_COMPILE__) || defined(__HIPCC__)
#include <hiprand/hiprand_kernel.h>
#endif
#define curandState                     hiprandState
#define curandState_t                   hiprandState_t
#define curand_init                     hiprand_init
#define curand_uniform                  hiprand_uniform
#define curand_normal                   hiprand_normal
#define curand_poisson                  hiprand_poisson

/*
 * Runtime verdict for fp32 hardware linear texture filtering (defined in
 * gpuUtils.cu). Cached and logged once per process.
 */
bool tigre_hw_linear_supported();

/*
 * Per-TU device flag mirroring the verdict, and the host helper that loads it.
 * Both are header-local (static), so each translation unit binds to its own
 * copy and no -fgpu-rdc is required: tex3D_TIGRE reads the flag in the same TU
 * whose host code called tigre_sync_hw_linear() after creating the texture.
 * Only emitted when a .cu is compiled by hipcc; plain-C++ TUs that include this
 * shim (e.g. GpuIds.cpp) never touch the texture path.
 */
#if defined(__HIPCC__)
static __constant__ int tigre_hw_linear_filter = 0;

/*
 * These are host helpers (called only from host code), but they must stay
 * declared in the device-compilation pass too: HIP single-source parses the
 * whole TU in both passes, so a host call site has to resolve the name even
 * though the host body is codegen'd only in the host pass.
 *
 * tex3D_TIGRE forwards to a single hardware tex3D<float> fetch when the flag is
 * set, and does the software 8-tap lerp when it is clear. So the flag is set
 * (1) for a Linear texture (hardware trilinear) AND for a deliberately Point
 * texture that the caller wants point-sampled as-is (matching CUDA, whose
 * tex3D_TIGRE is always the hardware fetch); it is cleared (0) only for the
 * Point fallback texture that needs software trilinear.
 */
__host__ static inline void tigre_set_hw_linear(int v) {
    hipMemcpyToSymbol(HIP_SYMBOL(tigre_hw_linear_filter), &v, sizeof(int));
}
__host__ static inline void tigre_sync_hw_linear() {
    tigre_set_hw_linear(tigre_hw_linear_supported() ? 1 : 0);
}
#endif /* __HIPCC__ */

#else /* CUDA */

#include <cuda_runtime_api.h>
#include <cuda.h>

/* On CUDA hardware linear filtering is always used; no self-test, no flag. */
static inline bool tigre_hw_linear_supported() { return true; }
static inline void tigre_sync_hw_linear() {}
static inline void tigre_set_hw_linear(int) {}

#endif /* USE_HIP */

/*
 * tex3D_TIGRE: replacement for the interpolated tex3D<float> reads.
 *
 * CUDA: forward to the hardware trilinear filter (the texture is created
 * cudaFilterModeLinear). HIP: when tigre_hw_linear_filter is set (the runtime
 * self-test confirmed hardware linear works and the texture was created
 * cudaFilterModeLinear) forward to the hardware fetch as well; otherwise the
 * texture is cudaFilterModePoint and we trilinearly interpolate in software.
 * CUDA's unnormalized texture
 * coordinate c samples between texels floor(c-0.5) and floor(c-0.5)+1 with
 * weight frac(c-0.5); a point sample of texel i is read at coordinate i+0.5.
 * Out-of-array neighbours read 0 through cudaAddressModeBorder, matching the
 * hardware filter's border behaviour.
 */
#if defined(__CUDACC__) || defined(__HIPCC__)
__device__ __forceinline__ float tex3D_TIGRE(cudaTextureObject_t tex, float x, float y, float z)
{
#if defined(USE_HIP)
    if (tigre_hw_linear_filter)
        return tex3D<float>(tex, x, y, z);

    const float xb = x - 0.5f, yb = y - 0.5f, zb = z - 0.5f;
    const float fx = floorf(xb), fy = floorf(yb), fz = floorf(zb);
    const float dx = xb - fx,   dy = yb - fy,   dz = zb - fz;

    const float c000 = tex3D<float>(tex, fx + 0.5f,        fy + 0.5f,        fz + 0.5f);
    const float c100 = tex3D<float>(tex, fx + 1.5f,        fy + 0.5f,        fz + 0.5f);
    const float c010 = tex3D<float>(tex, fx + 0.5f,        fy + 1.5f,        fz + 0.5f);
    const float c110 = tex3D<float>(tex, fx + 1.5f,        fy + 1.5f,        fz + 0.5f);
    const float c001 = tex3D<float>(tex, fx + 0.5f,        fy + 0.5f,        fz + 1.5f);
    const float c101 = tex3D<float>(tex, fx + 1.5f,        fy + 0.5f,        fz + 1.5f);
    const float c011 = tex3D<float>(tex, fx + 0.5f,        fy + 1.5f,        fz + 1.5f);
    const float c111 = tex3D<float>(tex, fx + 1.5f,        fy + 1.5f,        fz + 1.5f);

    const float c00 = c000 * (1.0f - dx) + c100 * dx;
    const float c10 = c010 * (1.0f - dx) + c110 * dx;
    const float c01 = c001 * (1.0f - dx) + c101 * dx;
    const float c11 = c011 * (1.0f - dx) + c111 * dx;

    const float c0 = c00 * (1.0f - dy) + c10 * dy;
    const float c1 = c01 * (1.0f - dy) + c11 * dy;

    return c0 * (1.0f - dz) + c1 * dz;
#else
    return tex3D<float>(tex, x, y, z);
#endif
}
#endif

#endif /* CUDA_TO_HIP_H */
