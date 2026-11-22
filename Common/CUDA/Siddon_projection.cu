/*-------------------------------------------------------------------------
 *
 * CUDA functions for ray-voxel intersection based projection
 *
 * This file has the necessary fucntiosn to perform X-ray CBCT projection
 * operation given a geaometry, angles and image. It usesthe so-called
 * Jacobs algorithm to compute efficiently the length of the x-rays over
 * voxel space.
 *
 * CODE by       Ander Biguri
 *               Sepideh Hatamikia (arbitrary rotation)
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
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "Siddon_projection.hpp"
#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include "launch_shapes.hpp"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                tigreCudaError("Ax:Siddon_projection", msg, cudaGetErrorString(__err));\
        } \
} while (0)
    
    
#define MAXTREADS 1024
#define PROJ_PER_BLOCK 9
#define MAX_PROJ_PER_LAUNCH (PROJ_PER_BLOCK*32)  // angles in one launch on small detectors (fan beam)
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
    
    void CreateTexture(const GpuIds& gpuids,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,bool alloc, cudaMemcpyKind kind);

__constant__ Point3D projParamsArrayDev[4*MAX_PROJ_PER_LAUNCH];  // Dev means it is on device


__global__ void vecAddInPlace(float *a, float *b, size_t n)
{
    size_t idx = (size_t)blockIdx.x*blockDim.x+threadIdx.x;
    // Make sure we do not go out of bounds
    if (idx < n)
        a[idx] = a[idx] + b[idx];
}

__global__ void kernelPixelDetector( Geometry geo,
        float* detector,
        const int currProjSetNumber,
        const int totalNoOfProjections,
        const int projPerLaunch,
        cudaTextureObject_t tex){
    
    
    unsigned long long u = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long v = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long long projNumber=blockIdx.z*blockDim.z+threadIdx.z;  // blocks of PROJ_PER_BLOCK angles stacked in z
    
    
    if (u>= geo.nDetecU || v>= geo.nDetecV || projNumber>=(unsigned long long)projPerLaunch)
        return;
    
#if IS_FOR_MATLAB_TIGRE
    size_t idx =  (size_t)(u * (unsigned long long)geo.nDetecV + v)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#else
    size_t idx =  (size_t)(v * (unsigned long long)geo.nDetecU + u)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#endif
    unsigned long indAlpha = currProjSetNumber*projPerLaunch+projNumber;  // This is the ABSOLUTE projection number in the projection array (for a given GPU)

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
    ///////
    // Siddon's ray-voxel intersection, optimized as in doi=10.1.1.55.7516
    //////
    // Also called Jacobs algorithms
    Point3D ray;
    // vector of Xray
    ray.x=pixel1D.x-source.x;
    ray.y=pixel1D.y-source.y;
    ray.z=pixel1D.z-source.z;
    float eps=0.001;
    ray.x=(fabsf(ray.x)<eps)? 0 : ray.x;
    ray.y=(fabsf(ray.y)<eps)? 0 : ray.y; 
    ray.z=(fabsf(ray.z)<eps)? 0 : ray.z; 
    // This variables are ommited because
    // bx,by,bz ={0,0,0}
    // dx,dy,dz ={1,1,1}
    // compute parameter values for x-ray parametric equation. eq(3-10)
    float axm,aym,azm;
    float axM,ayM,azM;
    // In the paper Nx= number of X planes-> Nvoxel+1
   
    axm=fminf(__fdividef(-source.x,ray.x),__fdividef(geo.nVoxelX-source.x,ray.x));
    aym=fminf(__fdividef(-source.y,ray.y),__fdividef(geo.nVoxelY-source.y,ray.y));
    azm=fminf(__fdividef(-source.z,ray.z),__fdividef(geo.nVoxelZ-source.z,ray.z));
    axM=fmaxf(__fdividef(-source.x,ray.x),__fdividef(geo.nVoxelX-source.x,ray.x));
    ayM=fmaxf(__fdividef(-source.y,ray.y),__fdividef(geo.nVoxelY-source.y,ray.y));
    azM=fmaxf(__fdividef(-source.z,ray.z),__fdividef(geo.nVoxelZ-source.z,ray.z));
    
    float am=fmaxf(fmaxf(axm,aym),azm);
    float aM=fminf(fminf(axM,ayM),azM);
    
    // line intersects voxel space ->   am<aM
    // The ray never enters the volume. RETURN HERE. Without it, execution
    // falls into the index block below, whose imin..kmax are only meaningful
    // for a ray that does enter. For a ray the eps guard has snapped to run
    // exactly along one voxel axis, with the source OUTSIDE the volume on
    // one of its zero axes, the two quotients for that axis share a sign, so
    // fminf and fmaxf collapse (axm == axM == -inf) and aM inherits -inf.
    // The other zero component then gives aM*ray = (-inf)*0 = NaN, Np is
    // NaN, and (unsigned long)NaN is 2^63 on CUDA: that thread iterates
    // ~9e18 times and the kernel never returns -- 100% SM, idle memory
    // controller, no error, no result. Measured: 2 rays in 261 million on an
    // odd x odd detector with an off-axis object; see
    // Python/tests/test_siddon_missing_ray.py and
    // Frontispiece/siddon_missed_ray_hang.png. A NEGATIVE Np does not do
    // this -- CUDA saturates negatives and -inf to 0, a zero-trip loop.
    if (am>=aM){
        detector[idx]=0;
        return;
    }

    // Walk the voxels the ray crosses, from where it enters the box (am) to where
    // it leaves (aM), adding up length x value.
    //
    // The first voxel is the one the ray is in just after entering. When the
    // entry point lies on a voxel boundary of another axis too, the direction
    // decides the side: moving up, floor(); moving down, ceil()-1. The next
    // crossing on each axis is then computed from that voxel's own faces, so the
    // two always agree. Before, the first voxel came from a midpoint and the
    // first crossings from ceil/floor of the entry point, separately; on such
    // rays they disagreed, giving a wrong first voxel or the whole ray one voxel
    // off along an axis. Rare for the whole image, but common when the image is
    // split in pieces, where many rays enter a piece through its z faces (#859).
    float pX=source.x+am*ray.x;
    float pY=source.y+am*ray.y;
    float pZ=source.z+am*ray.z;
    // 32-bit int: converting a 64-bit integer to float is very slow on the GPU
    // (64-bit indices made this kernel 5-8x slower, measured on an RTX 4070).
    int i=(ray.x<0.0f)? (int)ceilf(pX)-1 : (int)floorf(pX);
    int j=(ray.y<0.0f)? (int)ceilf(pY)-1 : (int)floorf(pY);
    int k=(ray.z<0.0f)? (int)ceilf(pZ)-1 : (int)floorf(pZ);
    i=(i<0)? 0 : ((i>(int)geo.nVoxelX-1)? (int)geo.nVoxelX-1 : i);
    j=(j<0)? 0 : ((j>(int)geo.nVoxelY-1)? (int)geo.nVoxelY-1 : j);
    k=(k<0)? 0 : ((k>(int)geo.nVoxelZ-1)? (int)geo.nVoxelZ-1 : k);

    // From here on, the current voxel is kept as the float coordinates of its
    // centre (index + 0.5), which is what the texture fetch takes: no
    // int-to-float conversion in the loop. Exact: integers + 0.5 are exact in
    // float up to 2^23 voxels per axis.
    float cX=(float)i+0.5f;
    float cY=(float)j+0.5f;
    float cZ=(float)k+0.5f;
    // eq(29), direction of update, and the offset from the centre to the face
    // ahead: +0.5 moving up, -0.5 moving down.
    float sX=(ray.x<0.0f)? -1.0f : 1.0f;
    float sY=(ray.y<0.0f)? -1.0f : 1.0f;
    float sZ=(ray.z<0.0f)? -1.0f : 1.0f;
    float hX=0.5f*sX;
    float hY=0.5f*sY;
    float hZ=0.5f*sZ;

    // Parameter of the next voxel boundary on each axis (infinite if the ray is
    // parallel to it): the boundary at coordinate b is at b*rX+oX. Each one is
    // computed from the voxel, not by adding a step to the previous one (eq 28):
    // in float32 those additions drift over the ~1000 boundaries of a large
    // image (errors of 2e-4 typical, 2e-3 at worst, at 512^3); this way the
    // error stays ~1e-5.
    float rX=(ray.x!=0.0f)? 1.0f/ray.x : 0.0f;
    float rY=(ray.y!=0.0f)? 1.0f/ray.y : 0.0f;
    float rZ=(ray.z!=0.0f)? 1.0f/ray.z : 0.0f;
    float oX=-source.x*rX;
    float oY=-source.y*rY;
    float oZ=-source.z*rZ;
    float ax=(ray.x!=0.0f)? fmaf(cX+hX,rX,oX) : INFINITY;
    float ay=(ray.y!=0.0f)? fmaf(cY+hY,rY,oY) : INFINITY;
    float az=(ray.z!=0.0f)? fmaf(cZ+hZ,rZ,oZ) : INFINITY;

    float maxlength=__fsqrt_rd(ray.x*ray.x*geo.dVoxelX*geo.dVoxelX+ray.y*ray.y*geo.dVoxelY*geo.dVoxelY+ray.z*ray.z*geo.dVoxelZ*geo.dVoxelZ);
    float sum=0.0f;
    float ac=am;
    // At most one step per voxel boundary, plus slack: the loop always ends, even
    // if rounding never quite lets a crossing reach aM.
    //
    // The texture fetch comes first and is the same for all threads; only the
    // update depends on which axis is crossed. (Master fetched inside each of
    // the three branches: when the threads of a warp crossed different axes,
    // the warp ran up to three fetches per step. This loop is ~1.5-1.8x faster.)
    const int maxsteps=(int)geo.nVoxelX+(int)geo.nVoxelY+(int)geo.nVoxelZ+3;
    for (int step=0; step<maxsteps; step++){
        float aminc=fminf(fminf(ax,ay),az);
        sum+=(fminf(aminc,aM)-ac)*tex3D<float>(tex, cX, cY, cZ);
        if (aminc>=aM)
            break;
        if (ax==aminc){
            cX+=sX; ac=ax; ax=fmaf(cX+hX,rX,oX);
        }else if(ay==aminc){
            cY+=sY; ac=ay; ay=fmaf(cY+hY,rY,oY);
        }else{
            cZ+=sZ; ac=az; az=fmaf(cZ+hZ,rZ,oZ);
        }
    }
    detector[idx]=sum*maxlength;
}


int siddon_ray_projection(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    cudaCheckErrors("Device query fail");
    // If a CUDA error raises (Python build) before the end of this function,
    // onError releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    if (gpuData && deviceCount != 1) {
        mexErrMsgIdAndTxt("Ax:GPUdata", "GPU arrays can only be projected on one GPU, the one they are on.");
    }
    // With GPU data the copies to and from the caller's array are device-to-device
    // and asynchronous, so they go on the compute stream: that orders them with
    // the kernels that reuse the double-buffered projections. (Host copies to
    // pageable memory are synchronous, so the copy stream is safe for them.)
    const int copyStream = gpuData ? 0 : 1;
    std::vector<int> devIds(deviceCount);  // copied: the release actions must not hold on to gpuids
    for (int d = 0; d < deviceCount; d++) { devIds[d] = gpuids[d]; onError.addDevice(devIds[d]); }
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection:GPUselect","There are no available device(s) that support CUDA\n");
    }
    //
    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning thrown)
    // Check the available devices, and if they are the same
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Ax:Siddon_projection:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;
    
    // Check free memory
    size_t mem_GPU_global;
    checkFreeMemory(gpuids, &mem_GPU_global);

    size_t mem_image=                 (unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
    size_t mem_proj=                  (unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV*sizeof(float);
    // Angles per kernel launch: PROJ_PER_BLOCK, or a multiple of it on small detectors (fan beam), whose launches of
    // PROJ_PER_BLOCK angles are too small to keep the GPU busy (tigre::projectionsPerLaunch).
    const int projPerLaunch=tigre::projectionsPerLaunch(geo.nDetecV, mem_proj, PROJ_PER_BLOCK, MAX_PROJ_PER_LAUNCH);
    
    // Does everything fit in the GPUs?
    bool fits_in_memory = mem_image+2*projPerLaunch*mem_proj<mem_GPU_global;
    unsigned int splits=1;
    if (!fits_in_memory) {
        // Nope nope.
        // approx free memory we have. We already have left some extra 5% free for internal stuff
        // we need a second projection memory to combine multi-GPU stuff.
        // The smallest working set is the 4 projection buffers and one image slice.
        const size_t mem_slice=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*sizeof(float);
        if (mem_GPU_global <= 4*projPerLaunch*mem_proj+mem_slice)
            tigreNotEnoughGpuMemory("Ax:Siddon_projection:GPUmemory","Forward projection (Siddon)",4*projPerLaunch*mem_proj+mem_slice,mem_GPU_global);
        size_t mem_free=mem_GPU_global-4*projPerLaunch*mem_proj;
        splits=mem_image/mem_free+1;// Ceil of the truncation
        if (splits>geo.nVoxelZ) splits=geo.nVoxelZ;  // at least one slice per piece
    }
    size_t num_bytes_proj = (size_t)projPerLaunch*geo.nDetecU*geo.nDetecV * sizeof(float);
    // The image texture (a 3D array on each GPU, for the largest piece) is allocated first. A 3D array needs more
    // GPU memory than its raw size, so it can fail where the estimate above said it fits: then split the image in
    // more pieces and try again, down to one slice per piece. The projection buffers allocated below are tried
    // with it (allocated and freed), so that they fit too.
    cudaArray **d_cuArrTex = new cudaArray*[deviceCount]();
    Geometry* geoArray = NULL;
    while (true){
        geoArray = (Geometry*)malloc(splits*sizeof(Geometry));
        splitImage(splits,geo,geoArray,nangles);
        const cudaExtent extent = make_cudaExtent(geoArray[0].nVoxelX, geoArray[0].nVoxelY, geoArray[0].nVoxelZ);
        if (tigre::tryMalloc3DArrays(devIds.data(), deviceCount, d_cuArrTex, extent, (splits>1 ? 4 : 2), num_bytes_proj, "Ax:Siddon_projection:GPU"))
            break;
        const unsigned long piece=geoArray[0].nVoxelZ;
        freeGeoArray(splits, geoArray);
        if (piece<=1 || splits>=geo.nVoxelZ){
            delete[] d_cuArrTex;
            size_t free_bytes=0, total_bytes=0;
            cudaMemGetInfo(&free_bytes, &total_bytes);
            tigreNotEnoughGpuMemory("Ax:Siddon_projection:GPUmemory","Forward projection (Siddon)",4*num_bytes_proj+(size_t)geo.nVoxelX*geo.nVoxelY*sizeof(float),free_bytes);
        }
        while (splits<geo.nVoxelZ && (geo.nVoxelZ+splits-1)/splits>=piece) splits++;  // smaller pieces
    }
    fits_in_memory = (splits==1);
    onError.add([geoArray, splits]{ freeGeoArray(splits, geoArray); });
    onError.add([d_cuArrTex, deviceCount, devIds]{
        for (int d = 0; d < deviceCount; d++) { cudaSetDevice(devIds[d]); if (d_cuArrTex[d]) cudaFreeArray(d_cuArrTex[d]); }
        delete[] d_cuArrTex;
    });
    
    // Allocate axuiliary memory for projections on the GPU to accumulate partial results
    float ** dProjection_accum = 0;
    if (!fits_in_memory){
        dProjection_accum=(float**)malloc(2*deviceCount*sizeof(float*));
        memset(dProjection_accum, 0, 2*deviceCount*sizeof(float*));
        onError.add([dProjection_accum, deviceCount, devIds]{
            for (int s = 0; s < 2*deviceCount; s++) { cudaSetDevice(devIds[s/2]); cudaFree(dProjection_accum[s]); }
            free(dProjection_accum);
        });
        for (dev = 0; dev < deviceCount; dev++) {
            cudaSetDevice(gpuids[dev]);
            for (int i = 0; i < 2; ++i){
                tigre::checkMalloc(cudaMalloc((void**)&dProjection_accum[dev*2+i], num_bytes_proj), "Ax:Siddon_projection:GPUmemory", "Forward projection (Siddon)", num_bytes_proj);
                cudaMemset(dProjection_accum[dev*2+i],0,num_bytes_proj); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
                cudaCheckErrors("cudaMallocauxiliarty projections fail");
            }
        }
    }
    
    // This is happening regarthless if the image fits on memory
    float** dProjection=(float**)malloc(2*deviceCount*sizeof(float*));
    memset(dProjection, 0, 2*deviceCount*sizeof(float*));
    onError.add([dProjection, deviceCount, devIds]{
        for (int s = 0; s < 2*deviceCount; s++) { cudaSetDevice(devIds[s/2]); cudaFree(dProjection[s]); }
        free(dProjection);
    });
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        
        for (int i = 0; i < 2; ++i){
            tigre::checkMalloc(cudaMalloc((void**)&dProjection[dev*2+i], num_bytes_proj), "Ax:Siddon_projection:GPUmemory", "Forward projection (Siddon)", num_bytes_proj);
            cudaMemset(dProjection[dev*2+i]  ,0,num_bytes_proj); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
            cudaCheckErrors("cudaMalloc projections fail");
        }
    }
    
    
    //Pagelock memory for synchronous copy.
    // Lets try to make the host memory pinned:
    // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
    if (gpuData) isHostRegisterSupported = 0;  // GPU memory: nothing to pin
    // empirical testing shows that when the image split is smaller than 1 (also implies the image is not very big), the time to
    // pin the memory is greater than the lost time in Synchronously launching the memcpys. This is only worth it when the image is too big.
#ifndef NO_PINNED_MEMORY
    // tigre::hostRegister: an array that is already page-locked (e.g. pinned by PyTorch) is used as it is and not unregistered here.
    bool img_registered=false;
    if (tigre::worthPageLocking(isHostRegisterSupported, gpuData, splits, deviceCount)){
        img_registered=tigre::hostRegister(img, (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*(size_t)sizeof(float));
        // The caller's array: left registered, the next call with it would fail.
        if (img_registered) onError.add([img]{ cudaHostUnregister(img); });
    }
#endif
    cudaCheckErrors("Error pinning memory");

    
    
    // auxiliary variables
    Point3D source, deltaU, deltaV, uvOrigin;
    Point3D* projParamsArrayHost = 0;
    cudaMallocHost((void**)&projParamsArrayHost,4*projPerLaunch*sizeof(Point3D));
    onError.add([projParamsArrayHost]{ cudaFreeHost(projParamsArrayHost); });
    cudaCheckErrors("Error allocating auxiliary constant memory");
    
    // Create Streams for overlapping memcopy and compute
    int nStreams=deviceCount*2;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    memset(stream, 0, nStreams*sizeof(cudaStream_t));
    onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
    
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        for (int i = 0; i < 2; ++i){
            cudaStreamCreate(&stream[i+dev*2]);
            
        }
    }
    cudaCheckErrors("Stream creation fail");

    int nangles_device=(nangles+deviceCount-1)/deviceCount;
    int nangles_last_device=(nangles-(deviceCount-1)*nangles_device);
    unsigned int noOfKernelCalls = (nangles_device+projPerLaunch-1)/projPerLaunch;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by projPerLaunch
    unsigned int noOfKernelCallsLastDev = (nangles_last_device+projPerLaunch-1)/projPerLaunch; // we will use this in the memory management.
    int projection_this_block;
    // Value-initialised (zeroed) so CreateTexture can tell an empty slot from one
    // holding a texture object it must destroy first.
    cudaTextureObject_t *texImg = new cudaTextureObject_t[deviceCount]();
    onError.add([texImg, deviceCount, devIds]{
        for (int d = 0; d < deviceCount; d++) {
            cudaSetDevice(devIds[d]);
            if (texImg[d]) cudaDestroyTextureObject(texImg[d]);
        }
        delete[] texImg;
    });
    
    for (unsigned int sp=0;sp<splits;sp++){
        
        // Create texture objects for all GPUs
        
        
        size_t linear_idx_start;
        //First one should always be  the same size as all the rest but the last
        linear_idx_start= (size_t)sp*(size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        
        
        CreateTexture(gpuids,&img[linear_idx_start],geoArray[sp],d_cuArrTex,texImg,false,(gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice));  // arrays allocated above
        cudaCheckErrors("Texture object creation fail");
        
        
        // Prepare kernel lauch variables
        
        int divU,divV;
        tigre::detectorBlockShape(geo.nDetecU, geo.nDetecV, PIXEL_SIZE_BLOCK, &divU, &divV);  // thin detectors (fan beam)
        dim3 grid((geoArray[sp].nDetecU+divU-1)/divU,(geoArray[0].nDetecV+divV-1)/divV,projPerLaunch/PROJ_PER_BLOCK);
        dim3 block(divU,divV,PROJ_PER_BLOCK);  // PROJ_PER_BLOCK angles per block, projPerLaunch/PROJ_PER_BLOCK blocks in z
        
        unsigned int proj_global;
        // Now that we have prepared the image (piece of image) and parameters for kernels
        // we project for all angles.
        for (unsigned int i=0; i<noOfKernelCalls; i++) {
            for (dev=0;dev<deviceCount;dev++){
                cudaSetDevice(gpuids[dev]);
                
                for(unsigned int j=0; j<projPerLaunch; j++){
                    proj_global=(i*projPerLaunch+j)+dev*nangles_device;
                    if (proj_global>=nangles)
                        break;
                    if ((i*projPerLaunch+j)>=nangles_device)
                        break;
                    geoArray[sp].alpha=angles[proj_global*3];
                    geoArray[sp].theta=angles[proj_global*3+1];
                    geoArray[sp].psi  =angles[proj_global*3+2];
                    
                    
                    //precomute distances for faster execution
                    //Precompute per angle constant stuff for speed
                    computeDeltas_Siddon(geoArray[sp],proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                    //Ray tracing!
                    projParamsArrayHost[4*j]=uvOrigin;		// 6*j because we have 6 Point3D values per projection
                    projParamsArrayHost[4*j+1]=deltaU;
                    projParamsArrayHost[4*j+2]=deltaV;
                    projParamsArrayHost[4*j+3]=source;
                    
                }
                cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*projPerLaunch,0,cudaMemcpyHostToDevice,stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2]);
                cudaCheckErrors("kernel fail");
                kernelPixelDetector<<<grid,block,0,stream[dev*2]>>>(geoArray[sp],dProjection[(i%2)+dev*2],i,nangles_device,projPerLaunch,texImg[dev]);
            }


            // Now that the computation is happening, we need to either prepare the memory for
            // combining of the projections (splits>1) and start removing previous results.
            
            
            // If our image does not fit in memory then we need to make sure we accumulate previous results too.
            // This is done in 2 steps: 
            // 1)copy previous results back into GPU 
            // 2)accumulate with current results
            // The code to take them out is the same as when there are no splits needed
            if( !fits_in_memory&&sp>0)
            {
                // 1) grab previous results and put them in the auxiliary variable dProjection_accum
                for (dev = 0; dev < deviceCount; dev++)
                {
                    cudaSetDevice(gpuids[dev]);
                    //Global index of FIRST projection on this set on this GPU
                    proj_global=i*projPerLaunch+dev*nangles_device;
                    if(proj_global>=nangles) 
                        break;

                    // Unless its the last projection set, we have projPerLaunch angles. Otherwise...
                    if(i+1==noOfKernelCalls) //is it the last block?
                        projection_this_block=min(nangles_device-(noOfKernelCalls-1)*projPerLaunch, //the remaining angles that this GPU had to do (almost never projPerLaunch)
                                                  nangles-proj_global);                              //or whichever amount is left to finish all (this is for the last GPU)
                    else
                        projection_this_block=projPerLaunch;

                    cudaMemcpyAsync(dProjection_accum[(i%2)+dev*2], result[proj_global], (size_t)projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice),stream[dev*2+copyStream]);
                }
                //  2) take the results from current compute call and add it to the code in execution.
                for (dev = 0; dev < deviceCount; dev++)
                {
                    cudaSetDevice(gpuids[dev]);
                    //Global index of FIRST projection on this set on this GPU
                    proj_global=i*projPerLaunch+dev*nangles_device;
                    if(proj_global>=nangles) 
                        break;

                    // Unless its the last projection set, we have projPerLaunch angles. Otherwise...
                    if(i+1==noOfKernelCalls) //is it the last block?
                        projection_this_block=min(nangles_device-(noOfKernelCalls-1)*projPerLaunch, //the remaining angles that this GPU had to do (almost never projPerLaunch)
                                                  nangles-proj_global);                              //or whichever amount is left to finish all (this is for the last GPU)
                    else
                        projection_this_block=projPerLaunch;

                    cudaStreamSynchronize(stream[dev*2+1]); // wait until copy is finished
                    vecAddInPlace<<<(unsigned int)(((size_t)geo.nDetecU*geo.nDetecV*projection_this_block+MAXTREADS-1)/MAXTREADS),MAXTREADS,0,stream[dev*2]>>>(dProjection[(i%2)+dev*2],dProjection_accum[(i%2)+dev*2],(size_t)geo.nDetecU*geo.nDetecV*projection_this_block);
                }
            } // end accumulation case, where the image needs to be split 

            // Now, lets get out the projections from the previous execution of the kernels.
            if (i>0){
                for (dev = 0; dev < deviceCount; dev++)
                {
                    cudaSetDevice(gpuids[dev]);
                    //Global index of FIRST projection on previous set on this GPU
                    proj_global=(i-1)*projPerLaunch+dev*nangles_device;
                    if (dev+1==deviceCount) {    //is it the last device?
                        // projections assigned to this device is >=nangles_device-(deviceCount-1) and < nangles_device
                        if (i-1 < noOfKernelCallsLastDev) {
                            // The previous set(block) was not empty.
                            projection_this_block=min(projPerLaunch, nangles-proj_global);
                        }
                        else {
                            // The previous set was empty.
                            // This happens if deviceCount > projPerLaunch+1.
                            // e.g. projPerLaunch = 9, deviceCount = 11, nangles = 199.
                            // e.g. projPerLaunch = 1, deviceCount =  3, nangles =   7.
                            break;
                        }
                    }
                    else {
                        projection_this_block=projPerLaunch;
                    }
                    cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(i%2))+dev*2],  (size_t)projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost),stream[dev*2+copyStream]);
                }
            }
            // Make sure Computation on kernels has finished before we launch the next batch.
            // ...and the copy-out of the batch before it (copy stream), which reads the
            // buffer the next kernel writes. Today that copy goes to pageable host memory
            // and is already finished when cudaMemcpyAsync returns, but nothing else
            // guarantees it (a page-locked result would race), see #860.
            for (dev = 0; dev < deviceCount; dev++){
                cudaSetDevice(gpuids[dev]);
                cudaStreamSynchronize(stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2+1]);
            }
        }
        
        
         // We still have the last set of projections to get out of GPUs
        for (dev = 0; dev < deviceCount; dev++)
        {
            cudaSetDevice(gpuids[dev]);
            //Global index of FIRST projection on this set on this GPU
            proj_global=(noOfKernelCalls-1)*projPerLaunch+dev*nangles_device;
            if(proj_global>=nangles) 
                break;
            // How many projections are left here?
            projection_this_block=min(nangles_device-(noOfKernelCalls-1)*projPerLaunch, //the remaining angles that this GPU had to do (almost never projPerLaunch)
                                      nangles-proj_global);                              //or whichever amount is left to finish all (this is for the last GPU)

            cudaDeviceSynchronize(); //Not really necessary, but just in case, we los nothing. 
            cudaCheckErrors("Error at copying the last set of projections out (or in the previous copy)");
            cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(noOfKernelCalls%2))+dev*2], (size_t)projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost),stream[dev*2+copyStream]);
        }
        // Make sure everyone has done their bussiness before the next image split:
        cudaDeviceSynchronize();
    } // End image split loop.
    
    cudaCheckErrors("Main loop  fail");
    onError.dismiss();  // from here on the function releases everything itself
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDestroyTextureObject(texImg[dev]);
            cudaFreeArray(d_cuArrTex[dev]);
    }
    delete[] texImg; texImg = 0;
    delete[] d_cuArrTex; d_cuArrTex = 0;
    // Freeing Stage
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaFree(dProjection[dev*2]);
        cudaFree(dProjection[dev*2+1]);
        
    }
    free(dProjection);
    
    if(!fits_in_memory){
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaFree(dProjection_accum[dev*2]);
            cudaFree(dProjection_accum[dev*2+1]);
            
        }
        free(dProjection_accum);
    }
    freeGeoArray(splits,geoArray);
    cudaFreeHost(projParamsArrayHost);
   
    
    for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]) ;
#ifndef NO_PINNED_MEMORY
    if (img_registered){
        cudaHostUnregister(img);
    }
    cudaCheckErrors("cudaFree  fail");
#endif
    //cudaDeviceReset();
    return 0;
}




void CreateTexture(const GpuIds& gpuids,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,bool alloc, cudaMemcpyKind kind)
{
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    const unsigned int num_devices = gpuids.GetLength();
    if(alloc){
        for (unsigned int dev = 0; dev < num_devices; dev++){
            cudaSetDevice(gpuids[dev]);
            
            //cudaArray Descriptor
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            //cuda Array
            tigre::checkMalloc(cudaMalloc3DArray(&d_cuArrTex[dev], &channelDesc, extent), "Ax:Siddon_projection:GPUmemory", "Forward projection (Siddon)", (size_t)extent.width*extent.height*extent.depth*sizeof(float));
        }
    }
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemcpy3DParms copyParams = {0};
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)imagedata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[dev];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
        cudaMemcpy3DAsync(&copyParams);
    }
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaSetDevice(gpuids[dev]);
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
        // Each image split creates a new texture object in this slot. Destroy the one
        // it replaces, otherwise every call that splits the image leaks one per split
        // (same failure as #617). The kernels that used it have finished: its array
        // has just been overwritten above.
        if (texImage[dev]) cudaDestroyTextureObject(texImage[dev]);
        cudaCreateTextureObject(&texImage[dev], &texRes, &texDescr, NULL);
        
    }
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaDeviceSynchronize();
    }
    cudaCheckErrors("Texture object creation fail");
}

/* This code generates the geometries needed to split the image properly in
 * cases where the entire image does not fit in the memory of the GPU
 **/
void splitImage(unsigned int splits,Geometry geo,Geometry* geoArray, unsigned int nangles){
    
    unsigned long splitsize=(geo.nVoxelZ+splits-1)/splits;// ceil if not divisible
    for(unsigned int sp=0;sp<splits;sp++){
        geoArray[sp]=geo;
        // All of them are splitsize, but the last one, possible
        geoArray[sp].nVoxelZ=((sp+1)*splitsize<geo.nVoxelZ)?  splitsize:  geo.nVoxelZ-splitsize*sp;
        geoArray[sp].sVoxelZ= geoArray[sp].nVoxelZ* geoArray[sp].dVoxelZ;
        
        // We need to redefine the offsets, as now each subimage is not aligned in the origin.
        geoArray[sp].offOrigZ=(float *)malloc(nangles*sizeof(float));
        for (unsigned int i=0;i<nangles;i++){
            geoArray[sp].offOrigZ[i]=geo.offOrigZ[i]-geo.sVoxelZ/2+sp*geoArray[0].sVoxelZ+geoArray[sp].sVoxelZ/2;
        }
        
    }
    
}

/* This code precomputes The location of the source and the Delta U and delta V (in the warped space)
 * to compute the locations of the x-rays. While it seems verbose and overly-optimized,
 * it does saves about 30% of each of the kernel calls. Thats something!
 **/
void computeDeltas_Siddon(Geometry geo,int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){

    
    Point3D S;
    S.x=geo.DSO[i];
    S.y=0;
    S.z=0;
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD[i]-geo.DSO[i]);   P.y  = geo.dDetecU*(-((double)geo.nDetecU/2.0)+0.5);       P.z  = geo.dDetecV*(((double)geo.nDetecV/2.0)-0.5);
    Pu0.x=0;                          Pu0.y= geo.dDetecU;                                    Pu0.z= 0;
    Pv0.x=0;                          Pv0.y= 0;                                              Pv0.z= geo.dDetecV*(-1);

    // Geometric transformations:
    // Now we have the Real world (OXYZ) coordinates of the bottom corner and its two neighbours.
    // The objective is to get a position of the detector in a coordinate system where:
    // 1-units are voxel size (in each direction can be different)
    // 2-The image has the its first voxel at (0,0,0)
    // 3-The image never rotates

    // To do that, we need to compute the "deltas" the detector, or "by how much
    // (in new xyz) does the voxels change when and index is added". To do that
    // several geometric steps needs to be changed

    //1.Roll,pitch,jaw
    // The detector can have a small rotation.
    // according to
    //"A geometric calibration method for cone beam CT systems" Yang K1, Kwan AL, Miller DF, Boone JM. Med Phys. 2006 Jun;33(6):1695-706.
    // Only the Z rotation will have a big influence in the image quality when they are small.
    // Still all rotations are supported

    // To roll pitch jaw, the detector has to be in centered in OXYZ.
    // NB: do not apply offsets to Pu0 and Pv0: they are directions, and are invariant through translations
    P.x=0;

    // Roll pitch yaw
    rollPitchYaw(geo,i,&P);
    rollPitchYaw(geo,i,&Pu0);
    rollPitchYaw(geo,i,&Pv0);
    //Now let's translate the points where they should be:
    // NB: do not apply offsets to Pu0 and Pv0: they are directions, and are invariant through translations
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);

    //1: Offset detector


    //S doesnt need to chagne


    //3: Rotate (around z)!
    Point3D Pfinal, Pfinalu0, Pfinalv0;
    Pfinal.x  =P.x;
    Pfinal.y  =P.y  +geo.offDetecU[i]; Pfinal.z  =P.z  +geo.offDetecV[i];
    Pfinalu0 = Pu0;
    Pfinalv0 = Pv0;

    eulerZYZ(geo,&Pfinal);
    eulerZYZ(geo,&Pfinalu0);
    eulerZYZ(geo,&Pfinalv0);
    eulerZYZ(geo,&S);

    //2: Offset image (instead of offseting image, -offset everything else)
    // NB: do not apply offsets to Pfinalu0 and Pfinalv0: they are directions, and are invariant through translations

    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    S.x=S.x-geo.offOrigX[i];               S.y=S.y-geo.offOrigY[i];               S.z=S.z-geo.offOrigZ[i];

    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2;
    S.x      =S.x+geo.sVoxelX/2;          S.y      =S.y+geo.sVoxelY/2;              S.z      =S.z      +geo.sVoxelZ/2;

    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S.x      =S.x/geo.dVoxelX;          S.y      =S.y/geo.dVoxelY;            S.z      =S.z/geo.dVoxelZ;


    //mexPrintf("COR: %f \n",geo.COR[i]);
    //5. apply COR. Wherever everything was, now its offesetd by a bit
    // NB: do not apply offsets to Pfinalu0 and Pfinalv0: they are directions, and are invariant through translations
    double CORx, CORy;
    CORx=-geo.COR[i]*sin(geo.alpha)/geo.dVoxelX;
    CORy= geo.COR[i]*cos(geo.alpha)/geo.dVoxelY;
    Pfinal.x+=CORx;   Pfinal.y+=CORy;
    S.x+=CORx; S.y+=CORy;

    // return

    *uvorigin=Pfinal;

    *deltaU=Pfinalu0;
    *deltaV=Pfinalv0;
    
    *source=S;
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
void rollPitchYaw(Geometry geo,int i, Point3D* point){
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x=cos(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.x
            +(cos(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) - sin(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.y
            +(cos(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) + sin(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.z;
    
    point->y=sin(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.x
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) + cos(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.y
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) - cos(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.z;
    
    point->z=-sin(geo.dPitch[i])*auxPoint.x
            +cos(geo.dPitch[i])*sin(geo.dYaw[i])*auxPoint.y
            +cos(geo.dPitch[i])*cos(geo.dYaw[i])*auxPoint.z;
    
}
void eulerZYZ(Geometry geo, Point3D* point){
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x=(+cos(geo.alpha)*cos(geo.theta)*cos(geo.psi)-sin(geo.alpha)*sin(geo.psi))*auxPoint.x+
            (-cos(geo.alpha)*cos(geo.theta)*sin(geo.psi)-sin(geo.alpha)*cos(geo.psi))*auxPoint.y+
            cos(geo.alpha)*sin(geo.theta)*auxPoint.z;
    
    point->y=(+sin(geo.alpha)*cos(geo.theta)*cos(geo.psi)+cos(geo.alpha)*sin(geo.psi))*auxPoint.x+
            (-sin(geo.alpha)*cos(geo.theta)*sin(geo.psi)+cos(geo.alpha)*cos(geo.psi))*auxPoint.y+
            sin(geo.alpha)*sin(geo.theta)*auxPoint.z;
    
    point->z=-sin(geo.theta)*cos(geo.psi)*auxPoint.x+
            sin(geo.theta)*sin(geo.psi)*auxPoint.y+
            cos(geo.theta)*auxPoint.z;
    
    
}
//______________________________________________________________________________
//
//      Function:       freeGeoArray
//
//      Description:    Frees the memory from the geometry array for multiGPU.
//______________________________________________________________________________
void freeGeoArray(unsigned int splits,Geometry* geoArray){
    for(unsigned int sp=0;sp<splits;sp++){
        free(geoArray[sp].offOrigZ);
    }
    free(geoArray);
}
//______________________________________________________________________________
//
//      Function:       checkFreeMemory
//
//      Description:    check available memory on devices
//______________________________________________________________________________
void checkFreeMemory(const GpuIds& gpuids, size_t *mem_GPU_global){
    size_t memfree;
    size_t memtotal;
    const int deviceCount = gpuids.GetLength();

    for (int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) *mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            mexWarnMsgIdAndTxt("Ax:Siddon_projection:GPUmemory","One (or more) of your GPUs is being heavily used by another program (less than half of its memory is free). TIGRE will try to run anyway, but it may need to split the problem into smaller pieces (slower) or run out of GPU memory.\n");
        }
        cudaCheckErrors("Check mem error");
        
        *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
    }
    *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
    
    //*mem_GPU_global= insert your known number here, in bytes.
}
#endif
