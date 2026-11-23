/*-------------------------------------------------------------------------
 *
 * CUDA functions for ray-voxel intersection based projection
 *
 * This file has the necessary fucntiosn to perform X-ray parallel projection
 * operation given a geaometry, angles and image. It usesthe so-called
 * Jacobs algorithm to compute efficiently the length of the x-rays over
 * voxel space. Its called Siddon because Jacobs algorithm its just a small
 * improvement over the traditional Siddons method.
 *
 * CODE by       Ander Biguri
 *
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


#include <algorithm>
#include <vector>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "Siddon_projection_parallel.hpp"
#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include "page_locking.hpp"
#include "launch_shapes.hpp"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                tigreCudaError("TIGRE:CUDA:Ax", msg, cudaGetErrorString(__err));\
        } \
} while (0)
    
    
// Declare the texture reference.
void CreateTextureParallel(const int* devs, int ndev, float* image, Geometry geo, cudaArray** d_cuArrTex, cudaTextureObject_t *texImage, cudaStream_t* stream, cudaMemcpyKind kind);


#define MAXTREADS 1024
#define PROJ_PER_BLOCK 9
#define PIXEL_SIZE_BLOCK 9
/*GEOMETRY DEFINITION
 *
 *                Detector plane, behind
 *            |-----------------------------|
 *            |                             |
 *            |                             |
 *            |                             |
 *            |                             |
 *            |      +--------+             |
 *            |     /        /|             |
 *   A Z      |    /        / |*D           |
 *   |        |   +--------+  |             |
 *   |        |   |        |  |             |
 *   |        |   |     *O |  +             |
 *    --->y   |   |        | /              |
 *  /         |   |        |/               |
 * V X        |   +--------+                |
 *            |-----------------------------|
 *
 *           *S
 *
 *
 *
 *
 *
 **/


__constant__ Point3D projParamsArrayDev[4*PROJ_PER_BLOCK];  // Dev means it is on device


__global__ void kernelPixelDetector_parallel( Geometry geo,
        float* detector, const int currProjSetNumber, const int totalNoOfProjections, cudaTextureObject_t tex){
    
    unsigned long long u = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long v = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long long projNumber=threadIdx.z;
            
    if (u>= geo.nDetecU || v>= geo.nDetecV || projNumber>=PROJ_PER_BLOCK)
        return;
    
    unsigned long indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    
#if IS_FOR_MATLAB_TIGRE
    size_t idx =  (size_t)(u  * (unsigned long long)geo.nDetecV + v)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#else
    size_t idx =  (size_t)(v  * (unsigned long long)geo.nDetecU + u)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#endif
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    

    /////// Get coordinates XYZ of pixel UV
    unsigned long pixelV = geo.nDetecV-v-1;
    unsigned long pixelU = u;
    Point3D pixel1D;
    pixel1D.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    pixel1D.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    pixel1D.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
    
    
    source.x=(source.x+pixelU*deltaU.x+pixelV*deltaV.x);
    source.y=(source.y+pixelU*deltaU.y+pixelV*deltaV.y);
    source.z=(source.z+pixelU*deltaU.z+pixelV*deltaV.z);
    ///////
    // Siddon's ray-voxel intersection, optimized as in doi=10.1.1.55.7516
    //////
    Point3D ray;
    // vector of Xray
    ray.x=pixel1D.x-source.x;
    ray.y=pixel1D.y-source.y;
    ray.z=pixel1D.z-source.z;

    float axm, aym, azm;
    float axM, ayM, azM;

    if (ray.x == 0.0f) {
        if (source.x >= 0.0f && source.x < geo.nVoxelX) {
            axm = -INFINITY;
            axM = INFINITY;
        } else {
            axm = INFINITY;
            axM = -INFINITY;
        }
    } else {
        float t1 = -source.x / ray.x;
        float t2 = (geo.nVoxelX - source.x) / ray.x;
        axm = fminf(t1, t2);
        axM = fmaxf(t1, t2);
    }

    if (ray.y == 0.0f) {
        if (source.y >= 0.0f && source.y < geo.nVoxelY) {
            aym = -INFINITY;
            ayM = INFINITY;
        } else {
            aym = INFINITY;
            ayM = -INFINITY;
        }
    } else {
        float t1 = -source.y / ray.y;
        float t2 = (geo.nVoxelY - source.y) / ray.y;
        aym = fminf(t1, t2);
        ayM = fmaxf(t1, t2);
    }

    if (ray.z == 0.0f) {
        if (source.z >= 0.0f && source.z < geo.nVoxelZ) {
            azm = -INFINITY;
            azM = INFINITY;
        } else {
            azm = INFINITY;
            azM = -INFINITY;
        }
    } else {
        float t1 = -source.z / ray.z;
        float t2 = (geo.nVoxelZ - source.z) / ray.z;
        azm = fminf(t1, t2);
        azM = fmaxf(t1, t2);
    }

    float am = fmaxf(fmaxf(axm, aym), azm);
    float aM = fminf(fminf(axM, ayM), azM);
    
    // line intersects voxel space ->   am<aM
    // See Siddon_projection.cu: without the return a missed ray falls
    // through to an Np built from indices that only mean something for a
    // hit; for an axis-parallel ray with the source outside the volume that
    // is (-inf)*0 = NaN, and (unsigned long)NaN = 2^63 spins forever. The
    // hang was NOT reproduced in this kernel (26 configurations tried); the
    // return here is by inspection.
    if (am>=aM){
        detector[idx]=0.0f;
        return;
    }

    // Compute max/min image INDEX for intersection eq(11-19)
    // Walk the voxels the ray crosses, from where it enters the box (am) to where
    // it leaves (aM), adding up length x value. Same traversal as the cone-beam
    // kernel (Siddon_projection.cu, #859):
    // - The first voxel is the one the ray is in just after entering; when the
    //   entry point lies on a voxel boundary, the direction decides the side.
    //   The next crossing on each axis is that voxel's own face ahead, so the
    //   two always agree.
    // - Each crossing is computed from its plane (Siddon's formula, one fmaf),
    //   not by adding a step to the previous one, which drifts in float32.
    // - The voxel is kept as the float coordinates of its centre, which the
    //   texture fetch takes directly; the fetch is the same for all threads,
    //   outside the branch on the axis crossed. Indices are 32-bit: 64-bit ones
    //   made this kernel several times slower.
    float pX=source.x+am*ray.x;
    float pY=source.y+am*ray.y;
    float pZ=source.z+am*ray.z;
    int i=(ray.x<0.0f)? (int)ceilf(pX)-1 : (int)floorf(pX);
    int j=(ray.y<0.0f)? (int)ceilf(pY)-1 : (int)floorf(pY);
    int k=(ray.z<0.0f)? (int)ceilf(pZ)-1 : (int)floorf(pZ);
    i=(i<0)? 0 : ((i>(int)geo.nVoxelX-1)? (int)geo.nVoxelX-1 : i);
    j=(j<0)? 0 : ((j>(int)geo.nVoxelY-1)? (int)geo.nVoxelY-1 : j);
    k=(k<0)? 0 : ((k>(int)geo.nVoxelZ-1)? (int)geo.nVoxelZ-1 : k);

    // Voxel centres (exact in float: integers + 0.5), direction, and the offset
    // from a centre to the face ahead.
    float cX=(float)i+0.5f;
    float cY=(float)j+0.5f;
    float cZ=(float)k+0.5f;
    float sX=(ray.x<0.0f)? -1.0f : 1.0f;
    float sY=(ray.y<0.0f)? -1.0f : 1.0f;
    float sZ=(ray.z<0.0f)? -1.0f : 1.0f;
    float hX=0.5f*sX;
    float hY=0.5f*sY;
    float hZ=0.5f*sZ;

    // Parameter of the next voxel boundary on each axis: the boundary at
    // coordinate b is at b*rX+oX (infinite if the ray is parallel to the axis).
    float rX=(ray.x!=0.0f)? 1.0f/ray.x : 0.0f;
    float rY=(ray.y!=0.0f)? 1.0f/ray.y : 0.0f;
    float rZ=(ray.z!=0.0f)? 1.0f/ray.z : 0.0f;
    float oX=-source.x*rX;
    float oY=-source.y*rY;
    float oZ=-source.z*rZ;
    float ax=(ray.x!=0.0f)? fmaf(cX+hX,rX,oX) : INFINITY;
    float ay=(ray.y!=0.0f)? fmaf(cY+hY,rY,oY) : INFINITY;
    float az=(ray.z!=0.0f)? fmaf(cZ+hZ,rZ,oZ) : INFINITY;

    float maxlength = sqrtf(ray.x * ray.x * geo.dVoxelX * geo.dVoxelX + ray.y * ray.y * geo.dVoxelY * geo.dVoxelY + ray.z * ray.z * geo.dVoxelZ * geo.dVoxelZ);
    float sum = 0.0f;
    float ac = am;
    // At most one step per voxel boundary, plus slack: the loop always ends, even
    // if rounding never quite lets a crossing reach aM.
    //
    // Four steps per pass, without a branch inside: each step adds its voxel, then
    // moves to the next voxel on the axis whose boundary comes first (selects, not
    // if/else). The fetch of each step does not wait for the sum of the one before,
    // so the four fetches are in flight together; before, every step waited for its
    // fetch before the next one could start, and a warp whose rays crossed different
    // axes ran the three branches one after the other. The sums are the same: after
    // the ray leaves (aminc>=aM), the remaining steps of the pass fetch a voxel each
    // but add nothing.
    const int maxsteps=(int)geo.nVoxelX+(int)geo.nVoxelY+(int)geo.nVoxelZ+3;
    bool out=false;
    for (int step=0; step<maxsteps; step+=4){
#pragma unroll
        for (int s=0; s<4; s++){
            const float aminc=fminf(fminf(ax,ay),az);
            const float t=tex3D<float>(tex, cX, cY, cZ);
            if (!out) sum+=(fminf(aminc,aM)-ac)*t;
            out=out || (aminc>=aM);
            const bool bx=(ax==aminc);
            const bool by=!bx && (ay==aminc);
            const bool bz=!bx && !by;
            ac=aminc;
            cX=bx? cX+sX : cX;  ax=bx? fmaf(cX+hX,rX,oX) : ax;
            cY=by? cY+sY : cY;  ay=by? fmaf(cY+hY,rY,oY) : ay;
            cZ=bz? cZ+sZ : cZ;  az=bz? fmaf(cZ+hZ,rZ,oZ) : az;
        }
        if (out)
            break;
    }
    detector[idx]=maxlength*sum;
}


// Adds b to a, element by element (the partial projections of the image pieces).
static __global__ void vecAddInPlaceParallel(float *a, const float *b, size_t n)
{
    size_t idx = (size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if (idx < n)
        a[idx] = a[idx] + b[idx];
}

// Free memory on the most used of these GPUs, less 5% for CUDA's own needs.
static size_t freeMemoryParallel(const int* devs, int ndev){
    size_t least=0;
    for (int d = 0; d < ndev; d++){
        size_t memfree=0, memtotal=0;
        cudaSetDevice(devs[d]);
        cudaMemGetInfo(&memfree,&memtotal);
        if (d==0 || memfree<least) least=memfree;
    }
    cudaCheckErrors("Check mem error");
    return (size_t)((double)least*0.95);
}

int siddon_ray_projection_parallel(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // Same scheme as the cone-beam projector (Siddon_projection.cu): the angles are shared between the GPUs, each GPU
    // holds the whole image (or, if it does not fit, one z-piece of it at a time, adding up the projections of the
    // pieces), and two projection buffers per GPU let the copy-out of one batch overlap the next (#286).
    int deviceCount = gpuids.GetLength();
    // If a CUDA error raises (Python build) before the end of this function,
    // onError releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection_parallel:GPUselect","There are no available device(s) that support CUDA\n");
    }
    if (gpuData && deviceCount != 1) {
        mexErrMsgIdAndTxt("Ax:GPUdata", "GPU arrays can only be projected on one GPU, the one they are on.");
    }
    // With GPU data the copies to and from the caller's array are device-to-device
    // and asynchronous, so they go on the compute stream: that orders them with
    // the kernels that reuse the double-buffered projections. (Host copies to
    // pageable memory are synchronous, so the copy stream is safe for them.)
    const int copyStream = gpuData ? 0 : 1;
    const cudaMemcpyKind kindIn  = gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice;
    const cudaMemcpyKind kindOut = gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost;
    std::vector<int> devIds(deviceCount);  // copied: the release actions must not hold on to gpuids
    for (int d = 0; d < deviceCount; d++) { devIds[d] = gpuids[d]; onError.addDevice(devIds[d]); }
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Ax:Siddon_projection_parallel:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;

    // Does the image fit on each GPU, with the projection buffers? If not, split it in z.
    const size_t mem_GPU_global=freeMemoryParallel(devIds.data(), deviceCount);
    const size_t mem_image=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*sizeof(float);
    const size_t num_bytes_proj=(size_t)PROJ_PER_BLOCK*(size_t)geo.nDetecU*(size_t)geo.nDetecV*sizeof(float);
    unsigned int splits=1;
    if (mem_image+2*num_bytes_proj>=mem_GPU_global){
        // Split: 2 more projection buffers to add up the pieces.
        const size_t mem_slice=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*sizeof(float);
        if (mem_GPU_global <= 4*num_bytes_proj+mem_slice)
            tigreNotEnoughGpuMemory("Ax:Siddon_projection_parallel:GPUmemory","Forward projection (Siddon, parallel beam)",4*num_bytes_proj+mem_slice,mem_GPU_global);
        splits=mem_image/(mem_GPU_global-4*num_bytes_proj)+1;
        if (splits>geo.nVoxelZ) splits=geo.nVoxelZ;  // at least one slice per piece
    }
    // The image texture (a 3D array on each GPU, for the largest piece) is allocated first, with the projection
    // buffers on top. A 3D array needs more GPU memory than its raw size, so it can fail where the estimate
    // above said it fits: then split the image in more pieces and try again, down to one slice per piece.
    // The projection buffers: two per GPU, and two more to add up the projections of the image pieces when the image
    // is split. tryMalloc3DArrays allocates them with the texture, each as its own allocation, so that the check that
    // everything fits covers them as they are used; they are handed out below.
    float** dProjection=(float**)malloc(2*deviceCount*sizeof(float*));
    memset(dProjection, 0, 2*deviceCount*sizeof(float*));
    float** dProjection_accum=(float**)malloc(2*deviceCount*sizeof(float*));
    memset(dProjection_accum, 0, 2*deviceCount*sizeof(float*));
    onError.add([dProjection, dProjection_accum, deviceCount, devIds]{
        for (int s = 0; s < 2*deviceCount; s++) { cudaSetDevice(devIds[s/2]); cudaFree(dProjection[s]); cudaFree(dProjection_accum[s]); }
        free(dProjection); free(dProjection_accum);
    });
    // GPU data on one GPU, image not split: the kernels write the projections straight into the caller's array, as
    // the backprojectors do, instead of into buffers that are then copied there. Only if it is one contiguous block
    // (it is, from Python), checked here. With a split image the pieces are added up in the buffers.
    bool outContiguous = gpuData && deviceCount == 1;
    for (int k = 1; outContiguous && k < nangles; k++)
        if (result[k] != result[0] + (size_t)k*geo.nDetecU*geo.nDetecV) outContiguous = false;
    std::vector<void*> projBuffers(4*(size_t)deviceCount, (void*)NULL);
    cudaArray **d_cuArrTex = new cudaArray*[deviceCount]();
    Geometry* geoArray = NULL;
    while (true){
        // No empty pieces: as many as the piece size needs (5 slices in 4 pieces of 2 would leave the last one empty).
        splits=(geo.nVoxelZ+(geo.nVoxelZ+splits-1)/splits-1)/((geo.nVoxelZ+splits-1)/splits);
        geoArray = (Geometry*)malloc(splits*sizeof(Geometry));
        splitImageInterp(splits,geo,geoArray,nangles);
        const cudaExtent extent = make_cudaExtent(geoArray[0].nVoxelX, geoArray[0].nVoxelY, geoArray[0].nVoxelZ);
        if (tigre::tryMalloc3DArrays(devIds.data(), deviceCount, d_cuArrTex, extent, (splits>1 ? 4 : (outContiguous ? 0 : 2)), num_bytes_proj, "Ax:Siddon_projection_parallel:GPU", projBuffers.data()))
            break;
        const unsigned long piece=geoArray[0].nVoxelZ;
        freeGeoArray(splits, geoArray);
        if (piece<=1 || splits>=geo.nVoxelZ){
            delete[] d_cuArrTex;
            size_t free_bytes=0, total_bytes=0;
            cudaMemGetInfo(&free_bytes, &total_bytes);
            tigreNotEnoughGpuMemory("Ax:Siddon_projection_parallel:GPUmemory","Forward projection (Siddon, parallel beam)",4*num_bytes_proj+(size_t)geo.nVoxelX*geo.nVoxelY*sizeof(float),free_bytes);
        }
        while (splits<geo.nVoxelZ && (geo.nVoxelZ+splits-1)/splits>=piece) splits++;  // smaller pieces
    }
    const bool fits_in_memory = (splits==1);
    onError.add([geoArray, splits]{ freeGeoArray(splits, geoArray); });
    onError.add([d_cuArrTex, deviceCount, devIds]{
        for (int d = 0; d < deviceCount; d++) { cudaSetDevice(devIds[d]); if (d_cuArrTex[d]) cudaFreeArray(d_cuArrTex[d]); }
        delete[] d_cuArrTex;
    });

    const bool directOut = outContiguous && fits_in_memory;
    {   // the projection buffers allocated with the texture (tryMalloc3DArrays)
        const int nb = fits_in_memory ? (directOut ? 0 : 2) : 4;
        for (dev = 0; dev < deviceCount && nb > 0; dev++){
            for (int i = 0; i < 2; ++i){
                dProjection[dev*2+i] = (float*)projBuffers[dev*nb+i];
                if (!fits_in_memory) dProjection_accum[dev*2+i] = (float*)projBuffers[dev*nb+2+i];
            }
        }
    }

    // Page-lock the image when it is copied in pieces or to several GPUs: the copies are then faster and asynchronous.
    bool registered=false;
#ifndef NO_PINNED_MEMORY
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,devIds[0]);
#endif
    if (tigre::worthPageLocking(isHostRegisterSupported, gpuData, splits, deviceCount)){
        registered=tigre::hostRegister(img, mem_image);
        // The caller's array: left registered, the next call with it would fail.
        if (registered) onError.add([img]{ cudaHostUnregister(img); });
    }
#endif
    cudaCheckErrors("Error pinning memory");

    Point3D source, deltaU, deltaV, uvOrigin;
    // The geometry constants of each launch are filled here and copied to __constant__ memory before the kernel.
    // Ordinary host memory: page-locking these few KB (allocating and freeing them with cudaMallocHost and
    // cudaFreeHost) cost ~0.45 ms per call and did not change the kernels' time (measured on an RTX 4070: the
    // kernels read the __constant__ copy either way).
    Point3D* projParamsArrayHost = 0;
    projParamsArrayHost=(Point3D*)malloc(4*PROJ_PER_BLOCK*sizeof(Point3D));
    if (!projParamsArrayHost) mexErrMsgIdAndTxt("TIGRE:hostMemory","Not enough host memory for the geometry constants");
    onError.add([projParamsArrayHost]{ free(projParamsArrayHost); });
    cudaCheckErrors("Error allocating auxiliary constant memory");

    // Two streams per GPU: compute (and, with GPU data, copies) on the first, host copies on the second.
    const int nStreams=deviceCount*2;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    memset(stream, 0, nStreams*sizeof(cudaStream_t));
    onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(devIds[dev]);
        for (int i = 0; i < 2; ++i)
            cudaStreamCreate(&stream[i+dev*2]);
    }
    cudaCheckErrors("Stream creation fail");
    // Clear the projection buffers on the stream of the kernels that write them (stream[dev*2]): that orders the
    // clearing before the kernels, so the host does not have to wait for it.
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(devIds[dev]);
        for (int i = 0; i < 2; ++i){
            if (dProjection[dev*2+i]) cudaMemsetAsync(dProjection[dev*2+i],0,num_bytes_proj,stream[dev*2]);
            if (!fits_in_memory) cudaMemsetAsync(dProjection_accum[dev*2+i],0,num_bytes_proj,stream[dev*2]);
        }
    }
    cudaCheckErrors("Projection buffer memset fail");

    // Value-initialised (zeroed) so CreateTextureParallel can tell an empty slot from one
    // holding a texture object it must destroy first.
    cudaTextureObject_t *texImg = new cudaTextureObject_t[deviceCount]();
    onError.add([texImg, deviceCount, devIds]{
        for (int d = 0; d < deviceCount; d++) {
            cudaSetDevice(devIds[d]);
            if (texImg[d]) cudaDestroyTextureObject(texImg[d]);
        }
        delete[] texImg;
    });

    const int nangles_device=(nangles+deviceCount-1)/deviceCount;
    const unsigned int noOfKernelCalls = (nangles_device+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK;  // the kernel checks the bounds if nangles_device is not divisible by PROJ_PER_BLOCK
    int divU,divV;
    tigre::detectorBlockShape(geo.nDetecU, geo.nDetecV, PIXEL_SIZE_BLOCK, &divU, &divV);  // thin detectors (2D)
    dim3 numBlocks((geo.nDetecU+divU-1)/divU,(geo.nDetecV+divV-1)/divV,1);
    dim3 threadsPerBlock(divU,divV,PROJ_PER_BLOCK);
    // Projections of batch i of this GPU (-1 if it has none): PROJ_PER_BLOCK, fewer for its last batch, fewer or
    // none on the last GPU.
    auto projectionsInBatch=[&](unsigned int i, int d)->int{
        const int first=(int)(i*PROJ_PER_BLOCK)+d*nangles_device;
        if (first>=nangles || (int)(i*PROJ_PER_BLOCK)>=nangles_device) return -1;
        int n=nangles_device-(int)(i*PROJ_PER_BLOCK);
        if (n>PROJ_PER_BLOCK) n=PROJ_PER_BLOCK;
        if (n>nangles-first) n=nangles-first;
        return n;
    };

    for (unsigned int sp=0;sp<splits;sp++){
        // The pieces are all the size of the first one but the last.
        const size_t linear_idx_start=(size_t)sp*(size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        CreateTextureParallel(devIds.data(),deviceCount,&img[linear_idx_start],geoArray[sp],d_cuArrTex,texImg,stream,kindIn);
        cudaCheckErrors("Texture object creation fail");

        for (unsigned int i=0; i<noOfKernelCalls; i++){
            for (dev=0;dev<deviceCount;dev++){
                cudaSetDevice(devIds[dev]);
                for(unsigned int j=0; j<PROJ_PER_BLOCK; j++){
                    const unsigned int proj_global=(i*PROJ_PER_BLOCK+j)+dev*nangles_device;
                    if (proj_global>=nangles || (i*PROJ_PER_BLOCK+j)>=nangles_device)
                        break;
                    geoArray[sp].alpha=angles[proj_global*3];
                    geoArray[sp].theta=angles[proj_global*3+1];
                    geoArray[sp].psi  =angles[proj_global*3+2];
                    //Precompute per angle constant stuff for speed
                    computeDeltas_Siddon_parallel(geoArray[sp],geoArray[sp].alpha,proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                    projParamsArrayHost[4*j]=uvOrigin;
                    projParamsArrayHost[4*j+1]=deltaU;
                    projParamsArrayHost[4*j+2]=deltaV;
                    projParamsArrayHost[4*j+3]=source;
                }
                cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2]);  // projParamsArrayHost is refilled for the next GPU
                kernelPixelDetector_parallel<<<numBlocks,threadsPerBlock,0,stream[dev*2]>>>(geoArray[sp],(directOut ? result[i*PROJ_PER_BLOCK] : dProjection[(i%2)+dev*2]),i,nangles_device,texImg[dev]);
            }
            cudaCheckErrors("Kernel fail");

            // Pieces after the first: add what the previous pieces gave for these projections (already in result).
            if (!fits_in_memory && sp>0){
                for (dev = 0; dev < deviceCount; dev++){
                    const int n=projectionsInBatch(i,dev);
                    if (n<0) break;
                    cudaSetDevice(devIds[dev]);
                    const unsigned int proj_global=i*PROJ_PER_BLOCK+dev*nangles_device;
                    cudaMemcpyAsync(dProjection_accum[(i%2)+dev*2], result[proj_global], (size_t)n*geo.nDetecV*geo.nDetecU*sizeof(float), kindIn, stream[dev*2+copyStream]);
                }
                for (dev = 0; dev < deviceCount; dev++){
                    const int n=projectionsInBatch(i,dev);
                    if (n<0) break;
                    cudaSetDevice(devIds[dev]);
                    cudaStreamSynchronize(stream[dev*2+1]);  // wait until the copy is finished
                    const size_t len=(size_t)geo.nDetecU*geo.nDetecV*n;
                    vecAddInPlaceParallel<<<(unsigned int)((len+MAXTREADS-1)/MAXTREADS),MAXTREADS,0,stream[dev*2]>>>(dProjection[(i%2)+dev*2],dProjection_accum[(i%2)+dev*2],len);
                }
            }

            // Copy out the previous batch while this one computes.
            if (i>0){
                for (dev = 0; dev < deviceCount; dev++){
                    const int n=projectionsInBatch(i-1,dev);
                    if (n<0) break;
                    cudaSetDevice(devIds[dev]);
                    const unsigned int proj_global=(i-1)*PROJ_PER_BLOCK+dev*nangles_device;
                    if (!directOut) cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(i%2))+dev*2], (size_t)n*geo.nDetecV*geo.nDetecU*sizeof(float), kindOut, stream[dev*2+copyStream]);
                }
            }
            // Wait for this batch's kernels, and for the copy-out of the batch before it (copy stream), which reads
            // the buffer the next kernel writes (#860).
            for (dev = 0; dev < deviceCount; dev++){
                cudaSetDevice(devIds[dev]);
                cudaStreamSynchronize(stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2+1]);
            }
        }

        // The last batch of each GPU.
        for (dev = 0; dev < deviceCount; dev++){
            const int n=projectionsInBatch(noOfKernelCalls-1,dev);
            if (n<0) break;
            cudaSetDevice(devIds[dev]);
            const unsigned int proj_global=(noOfKernelCalls-1)*PROJ_PER_BLOCK+dev*nangles_device;
            if (!directOut) cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(noOfKernelCalls%2))+dev*2], (size_t)n*geo.nDetecV*geo.nDetecU*sizeof(float), kindOut, stream[dev*2+copyStream]);
        }
        // Everything done before the next piece (it reads these results back) or before returning.
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(devIds[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("Projection fail");
    }

    onError.dismiss();  // from here on the function releases everything itself
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(devIds[dev]);
        cudaDestroyTextureObject(texImg[dev]);
        cudaFreeArray(d_cuArrTex[dev]);
        cudaFree(dProjection[dev*2]);
        cudaFree(dProjection[dev*2+1]);
        if (!fits_in_memory){
            cudaFree(dProjection_accum[dev*2]);
            cudaFree(dProjection_accum[dev*2+1]);
        }
        cudaStreamDestroy(stream[dev*2]);
        cudaStreamDestroy(stream[dev*2+1]);
    }
    delete[] texImg;
    delete[] d_cuArrTex;
    free(dProjection);
    free(dProjection_accum);
    free(stream);
    freeGeoArray(splits,geoArray);
    free(projParamsArrayHost);
    if (registered) cudaHostUnregister(img);
    cudaCheckErrors("cudaFree fail");
    return 0;
}



/* This code precomputes The location of the source and the Delta U and delta V (in the warped space)
 * to compute the locations of the x-rays. While it seems verbose and overly-optimized,
 * it does saves about 30% of each of the kernel calls. Thats something!
 **/
void computeDeltas_Siddon_parallel(Geometry geo, float angles,int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){
    Point3D S;
    
    S.x  =geo.DSO[i];   S.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       S.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD[i]-geo.DSO[i]);   P.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       P.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pu0.x=-(geo.DSD[i]-geo.DSO[i]);   Pu0.y= geo.dDetecU*(1-((float)geo.nDetecU/2)+0.5);       Pu0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pv0.x=-(geo.DSD[i]-geo.DSO[i]);   Pv0.y= geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       Pv0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-1);
    // Geometric trasnformations:
    P.x=0;Pu0.x=0;Pv0.x=0;
    
    // Roll pitch yaw
    rollPitchYaw(geo,i,&P);
    rollPitchYaw(geo,i,&Pu0);
    rollPitchYaw(geo,i,&Pv0);
    //Now lets translate the points where they should be:
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);
    Pu0.x=Pu0.x-(geo.DSD[i]-geo.DSO[i]);
    Pv0.x=Pv0.x-(geo.DSD[i]-geo.DSO[i]);

    S.x=0;
    // Roll pitch yaw
    rollPitchYaw(geo,i,&S);
    //Now lets translate the points where they should be:
    S.x=S.x+geo.DSO[i];

    //1: Offset detector
    
    //P.x
    P.y  =P.y  +geo.offDetecU[i];    P.z  =P.z  +geo.offDetecV[i];
    Pu0.y=Pu0.y+geo.offDetecU[i];    Pu0.z=Pu0.z+geo.offDetecV[i];
    Pv0.y=Pv0.y+geo.offDetecU[i];    Pv0.z=Pv0.z+geo.offDetecV[i];
    //S doesnt need to chagne
    
    
    //3: Rotate using full Euler angles (ZYZ convention: alpha, theta, psi)
    Point3D Pfinal = P;
    Point3D Pfinalu0 = Pu0;
    Point3D Pfinalv0 = Pv0;
    Point3D S2 = S;
    
    eulerZYZ(geo, &Pfinal);
    eulerZYZ(geo, &Pfinalu0);
    eulerZYZ(geo, &Pfinalv0);
    eulerZYZ(geo, &S2);
    
    //2: Offset image (instead of offseting image, -offset everything else)
    
    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    Pfinalu0.x=Pfinalu0.x-geo.offOrigX[i];   Pfinalu0.y=Pfinalu0.y-geo.offOrigY[i];   Pfinalu0.z=Pfinalu0.z-geo.offOrigZ[i];
    Pfinalv0.x=Pfinalv0.x-geo.offOrigX[i];   Pfinalv0.y=Pfinalv0.y-geo.offOrigY[i];   Pfinalv0.z=Pfinalv0.z-geo.offOrigZ[i];
    S2.x=S2.x-geo.offOrigX[i];               S2.y=S2.y-geo.offOrigY[i];               S2.z=S2.z-geo.offOrigZ[i];
    
    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2;
    Pfinalu0.x=Pfinalu0.x+geo.sVoxelX/2;    Pfinalu0.y=Pfinalu0.y+geo.sVoxelY/2;        Pfinalu0.z=Pfinalu0.z+geo.sVoxelZ/2;
    Pfinalv0.x=Pfinalv0.x+geo.sVoxelX/2;    Pfinalv0.y=Pfinalv0.y+geo.sVoxelY/2;        Pfinalv0.z=Pfinalv0.z+geo.sVoxelZ/2;
    S2.x      =S2.x+geo.sVoxelX/2;          S2.y      =S2.y+geo.sVoxelY/2;              S2.z      =S2.z      +geo.sVoxelZ/2;
    
    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S2.x      =S2.x/geo.dVoxelX;          S2.y      =S2.y/geo.dVoxelY;            S2.z      =S2.z/geo.dVoxelZ;
    
    
    
    //5. apply COR. Wherever everything was, now its offesetd by a bit
    float CORx, CORy;
    CORx=-geo.COR[i]*sin(geo.alpha)/geo.dVoxelX;
    CORy= geo.COR[i]*cos(geo.alpha)/geo.dVoxelY;
    Pfinal.x+=CORx;   Pfinal.y+=CORy;
    Pfinalu0.x+=CORx;   Pfinalu0.y+=CORy;
    Pfinalv0.x+=CORx;   Pfinalv0.y+=CORy;
    S2.x+=CORx; S2.y+=CORy;
    
    // return
    
    *uvorigin=Pfinal;
    
    deltaU->x=Pfinalu0.x-Pfinal.x;
    deltaU->y=Pfinalu0.y-Pfinal.y;
    deltaU->z=Pfinalu0.z-Pfinal.z;
    
    deltaV->x=Pfinalv0.x-Pfinal.x;
    deltaV->y=Pfinalv0.y-Pfinal.y;
    deltaV->z=Pfinalv0.z-Pfinal.z;
    
    *source=S2;
}
void CreateTextureParallel(const int* devs, int ndev, float* image, Geometry geo, cudaArray** d_cuArrTex, cudaTextureObject_t *texImage, cudaStream_t* stream, cudaMemcpyKind kind){
    // Copy this image piece into each GPU's 3D array and make a texture object of it. The arrays were allocated for
    // the largest piece; for a smaller (last) piece the array is allocated again at its size, so that the texture's
    // border, not slices left over from the previous piece, lies beyond its last slice.
    const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    for (int dev = 0; dev < ndev; dev++){
        cudaSetDevice(devs[dev]);
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        cudaExtent current = make_cudaExtent(0,0,0);
        if (d_cuArrTex[dev]) { cudaChannelFormatDesc d; unsigned int flags; cudaArrayGetInfo(&d, &current, &flags, d_cuArrTex[dev]); }
        if (current.depth != extent.depth || current.width != extent.width || current.height != extent.height){
            if (texImage[dev]) { cudaDestroyTextureObject(texImage[dev]); texImage[dev] = 0; }
            if (d_cuArrTex[dev]) { cudaFreeArray(d_cuArrTex[dev]); d_cuArrTex[dev] = 0; }
            tigre::checkMalloc(cudaMalloc3DArray(&d_cuArrTex[dev], &channelDesc, extent), "Ax:Siddon_projection_parallel:GPUmemory", "Forward projection (Siddon, parallel beam)", (size_t)extent.width*extent.height*extent.depth*sizeof(float));
        }
        cudaMemcpy3DParms copyParams = {0};
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)image, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[dev];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
        cudaMemcpy3DAsync(&copyParams,stream[dev*2]);  // the kernels' stream: ordered before them
    }
    for (int dev = 0; dev < ndev; dev++){
        cudaSetDevice(devs[dev]);
        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArrTex[dev];
        cudaTextureDesc     texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = false;
        texDescr.filterMode = cudaFilterModePoint;
        texDescr.addressMode[0] = cudaAddressModeBorder;
        texDescr.addressMode[1] = cudaAddressModeBorder;
        texDescr.addressMode[2] = cudaAddressModeBorder;
        texDescr.readMode = cudaReadModeElementType;
        // Each image piece creates a new texture object in this slot: destroy the one it replaces (#617). The
        // kernels that used it have finished.
        if (texImage[dev]) cudaDestroyTextureObject(texImage[dev]);
        cudaCreateTextureObject(&texImage[dev], &texRes, &texDescr, NULL);
    }
    cudaCheckErrors("Texture object creation fail");
}

#ifndef PROJECTION_HPP

float maxDistanceCubeXY(Geometry geo, float alpha,int i){
    ///////////
    // Compute initial "t" so we access safely as less as out of bounds as possible.
    //////////
    
    
    float maxCubX,maxCubY;
    // Forgetting Z, compute max distance: diagonal+offset
    maxCubX=(geo.sVoxelX/2+ abs(geo.offOrigX[i]))/geo.dVoxelX;
    maxCubY=(geo.sVoxelY/2+ abs(geo.offOrigY[i]))/geo.dVoxelY;
    
    return geo.DSO[i]/geo.dVoxelX-sqrt(maxCubX*maxCubX+maxCubY*maxCubY);
    
}

#endif
