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
 --------------------------------------------------------------------------
 Shared helpers for the SYCL backend: device list, one in-order queue per
 device, and the texture fetches the CUDA kernels get from texture hardware.
 */

#ifndef TIGRE_SYCL_HPP
#define TIGRE_SYCL_HPP

#include <sycl/sycl.hpp>
#include <vector>
#include <cstddef>

//! @brief SYCL devices TIGRE runs on, in device-id order.
//! GPUs are listed if there are any (Level Zero preferred over OpenCL, so a
//! GPU is not listed twice); otherwise CPU devices (is_cpu()), printing once
//! that TIGRE is running on the CPU. With TIGRE_SYCL_ALL_DEVICES=1, CPU devices
//! are listed after the GPUs. The usual ONEAPI_DEVICE_SELECTOR filter applies.
const std::vector<sycl::device>& tigreSyclDevices();

//! @brief In-order queue for device id iDeviceId (an index of tigreSyclDevices()).
sycl::queue& tigreSyclQueue(int iDeviceId);

//! @brief Largest single allocation in bytes: the device's max_mem_alloc_size,
//! optionally capped by TIGRE_SYCL_MAX_ALLOC_MB (for tests).
size_t tigreSyclMaxAlloc(int iDeviceId);

//! @brief Total memory budget in bytes: 95% of global memory, optionally
//! capped by TIGRE_SYCL_MEMORY_LIMIT_MB (for tests). This is not limited by
//! max_mem_alloc_size; image pieces are split so each buffer stays within
//! tigreSyclMaxAlloc().
size_t tigreSyclMemoryBudget(int iDeviceId);

//! @brief Free device memory in bytes, as far as the SYCL runtime reports it
//! (sycl::ext::intel::info::device::free_memory, e.g. on Level Zero GPUs), and
//! never more than tigreSyclMemoryBudget(). Devices that do not report it
//! (e.g. the OpenCL CPU device) give tigreSyclMemoryBudget().
size_t tigreSyclFreeMemory(int iDeviceId);

//! @brief Total device memory in bytes (global memory size).
size_t tigreSyclTotalMemory(int iDeviceId);

//! @brief tigreSyclFreeMemory() for TIGRE's memory checks (checkFreeMemory).
//! As in the CUDA build, warns (tag) if less than half of the device memory is
//! free, when the runtime reports free memory.
size_t tigreSyclCheckFreeMemory(int iDeviceId, const char* tag);

//! @brief How many image pieces are needed so the working set fits in
//! mem_budget and each single buffer is at most max_alloc. Returns 1 if no
//! split is required, or 0 if even one image slice plus the projection
//! buffers cannot fit (the caller should raise tigreNotEnoughGpuMemory).
unsigned int tigreSyclImageSplits(size_t mem_image, size_t mem_slice,
                                  size_t mem_proj_block, unsigned int nVoxelZ,
                                  size_t mem_budget, size_t max_alloc);

//! @brief If TIGRE_SYCL_SPLIT_LOG is set and is not "0", print the number of
//! image pieces used (TIGRE_SYCL_SPLITS=N) to stderr, for tests.
void tigreSyclLogSplits(unsigned int splits);

//! @brief Point-sampled fetch with border (zero) addressing, as a CUDA
//! cudaFilterModePoint + cudaAddressModeBorder texture with unnormalized
//! coordinates.
inline float tigreFetchPoint(const float* vol, long nx, long ny, long nz, long i, long j, long k){
    if (i<0 || j<0 || k<0 || i>=nx || j>=ny || k>=nz)
        return 0.0f;
    return vol[(size_t)k*(size_t)nx*(size_t)ny+(size_t)j*(size_t)nx+(size_t)i];
}

//! @brief Point-sampled fetch at unnormalized float coordinates (texel floor(x)).
inline float tigreFetchPoint3D(const float* vol, long nx, long ny, long nz, float x, float y, float z){
    return tigreFetchPoint(vol, nx, ny, nz, (long)sycl::floor(x), (long)sycl::floor(y), (long)sycl::floor(z));
}

//! @brief Trilinear fetch with border (zero) addressing at unnormalized
//! coordinates, as a CUDA cudaFilterModeLinear texture: texel centres sit at
//! integer+0.5.
inline float tigreFetchLinear3D(const float* vol, long nx, long ny, long nz, float x, float y, float z){
    float xb=x-0.5f, yb=y-0.5f, zb=z-0.5f;
    float fx=sycl::floor(xb), fy=sycl::floor(yb), fz=sycl::floor(zb);
    long i=(long)fx, j=(long)fy, k=(long)fz;
    float a=xb-fx, b=yb-fy, c=zb-fz;
    float v000=tigreFetchPoint(vol,nx,ny,nz,i  ,j  ,k  );
    float v100=tigreFetchPoint(vol,nx,ny,nz,i+1,j  ,k  );
    float v010=tigreFetchPoint(vol,nx,ny,nz,i  ,j+1,k  );
    float v110=tigreFetchPoint(vol,nx,ny,nz,i+1,j+1,k  );
    float v001=tigreFetchPoint(vol,nx,ny,nz,i  ,j  ,k+1);
    float v101=tigreFetchPoint(vol,nx,ny,nz,i+1,j  ,k+1);
    float v011=tigreFetchPoint(vol,nx,ny,nz,i  ,j+1,k+1);
    float v111=tigreFetchPoint(vol,nx,ny,nz,i+1,j+1,k+1);
    return (1-c)*((1-b)*((1-a)*v000+a*v100)+b*((1-a)*v010+a*v110))
          +   c *((1-b)*((1-a)*v001+a*v101)+b*((1-a)*v011+a*v111));
}

#endif // TIGRE_SYCL_HPP
