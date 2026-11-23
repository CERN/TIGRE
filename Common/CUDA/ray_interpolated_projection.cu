/*-------------------------------------------------------------------------
 *
 * CUDA functions for texture-memory interpolation based projection
 *
 * This file has the necessary fucntiosn to perform X-ray CBCT projection
 * operation given a geaometry, angles and image. It uses the 3D texture
 * memory linear interpolation to uniformily sample a path to integrate the
 * X-rays.
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
#include "ray_interpolated_projection.hpp"
#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include "page_locking.hpp"
#include "launch_shapes.hpp"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                TIGRE_MATLAB_ONLY(cudaDeviceReset());\
                tigreCudaError("TIGRE:Ax:interpolated", msg, cudaGetErrorString(__err));\
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
    void CreateTextureInterp(const GpuIds& gpuids,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,bool allocate, cudaMemcpyKind kind, cudaStream_t* stream);
__constant__ Point3D projParamsArrayDev[4*MAX_PROJ_PER_LAUNCH];  // Dev means it is on device
__constant__ float projFloatsArrayDev[2*MAX_PROJ_PER_LAUNCH];  // Dev means it is on device


__global__ void vecAddInPlaceInterp(float *a, float *b, size_t n)
{
    size_t idx = (size_t)blockIdx.x*blockDim.x+threadIdx.x;
    // Make sure we do not go out of bounds
    if (idx < n)
        a[idx] = a[idx] + b[idx];
}


template<bool sphericalrotation>
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

    unsigned long indAlpha = currProjSetNumber*projPerLaunch+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    
    float DSO = projFloatsArrayDev[2*projNumber+0];
    float cropdist_init = projFloatsArrayDev[2*projNumber+1];
    
    
    
    /////// Get coordinates XYZ of pixel UV
    unsigned long pixelV = geo.nDetecV-v-1;
    unsigned long pixelU = u;
    
    
    float vectX,vectY,vectZ;
    Point3D P;
    P.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    P.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    P.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
    
    // Length is the ray length in normalized space
    float length=__fsqrt_rd((source.x-P.x)*(source.x-P.x)+(source.y-P.y)*(source.y-P.y)+(source.z-P.z)*(source.z-P.z));
    //now legth is an integer of Nsamples that are required on this line
    length=ceilf(__fdividef(length,geo.accuracy));//Divide the directional vector by an integer
    vectX=__fdividef(P.x -source.x,length);
    vectY=__fdividef(P.y -source.y,length);
    vectZ=__fdividef(P.z -source.z,length);
    
    
//     //Integrate over the line
    float tx,ty,tz;
    float sum=0;
    float i;
    
    
    
//  Because I have no idea how to efficiently cutoff the legth path in 3D, a very upper limit is computed (see maxdistanceCuboid)
//  for the 3D case. However it would be bad to lose performance in the 3D case
//  TODO: can ge really improve this?
    if (sphericalrotation){
        if ((2*DSO/fminf(fminf(geo.dVoxelX,geo.dVoxelY),geo.dVoxelZ)+cropdist_init)/geo.accuracy  <   length)
            length=ceilf((2*DSO/fminf(fminf(geo.dVoxelX,geo.dVoxelY),geo.dVoxelZ)+cropdist_init)/geo.accuracy);
    }
    else{
        if ((2*DSO/fminf(geo.dVoxelX,geo.dVoxelY)+cropdist_init)/geo.accuracy  <   length)
            length=ceilf((2*DSO/fminf(geo.dVoxelX,geo.dVoxelY)+cropdist_init)/geo.accuracy);
    }
    
    
    //Length is not actually a length, but the amount of memreads with given accuracy ("samples per voxel")
    for (i=floorf(cropdist_init/geo.accuracy); i<=length; i=i+1){
        tx=vectX*i+source.x;
        ty=vectY*i+source.y;
        tz=vectZ*i+source.z;
        
        sum += tex3D<float>(tex, tx+0.5f, ty+0.5f, tz+0.5f); // this line is 94% of time.
    }
    
    float deltalength=sqrtf((vectX*geo.dVoxelX)*(vectX*geo.dVoxelX)+
            (vectY*geo.dVoxelY)*(vectY*geo.dVoxelY)+
            (vectZ*geo.dVoxelZ)*(vectZ*geo.dVoxelZ) );
    
    detector[idx]=sum*deltalength;
}



// legnth(angles)=3 x nagnles, as we have roll, pitch, yaw.
int interpolation_projection(float  *  img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    
    
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
        mexErrMsgIdAndTxt("Ax:Interpolated_projection:GPUselect","There are no available device(s) that support CUDA\n");
    }
    //
    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning thrown)
    // Check the available devices, and if they are the same
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Ax:Interpolated_projection:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;
    
    // Check free memory
    size_t mem_GPU_global;
    checkFreeMemory(gpuids,&mem_GPU_global);

    // printf("geo.nDetec (U, V) = %d, %d\n", geo.nDetecU, geo.nDetecV);
    
    size_t mem_image=(unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
    size_t mem_proj =(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV * sizeof(float);
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
            tigreNotEnoughGpuMemory("Ax:Interpolated_projection:GPUmemory","Forward projection (interpolated)",4*projPerLaunch*mem_proj+mem_slice,mem_GPU_global);
        size_t mem_free=mem_GPU_global-4*projPerLaunch*mem_proj;
        splits=mem_image/mem_free+1;// Ceil of the truncation
        if (splits>geo.nVoxelZ) splits=geo.nVoxelZ;  // at least one slice per piece
    }
    size_t num_bytes_proj = (size_t)projPerLaunch*geo.nDetecU*geo.nDetecV * sizeof(float);
    // The image texture (a 3D array on each GPU, for the largest piece) is allocated first. A 3D array needs more
    // GPU memory than its raw size, so it can fail where the estimate above said it fits: then split the image in
    // more pieces and try again, down to one slice per piece.
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
        if (tigre::tryMalloc3DArrays(devIds.data(), deviceCount, d_cuArrTex, extent, (splits>1 ? 4 : (outContiguous ? 0 : 2)), num_bytes_proj, "Ax:Interpolated_projection:GPU", projBuffers.data()))
            break;
        const unsigned long piece=geoArray[0].nVoxelZ;
        freeGeoArray(splits, geoArray);
        if (piece<=1 || splits>=geo.nVoxelZ){
            delete[] d_cuArrTex;
            size_t free_bytes=0, total_bytes=0;
            cudaMemGetInfo(&free_bytes, &total_bytes);
            tigreNotEnoughGpuMemory("Ax:Interpolated_projection:GPUmemory","Forward projection (interpolated)",4*num_bytes_proj+(size_t)geo.nVoxelX*geo.nVoxelY*sizeof(float),free_bytes);
        }
        while (splits<geo.nVoxelZ && (geo.nVoxelZ+splits-1)/splits>=piece) splits++;  // smaller pieces
    }
    fits_in_memory = (splits==1);
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
    if (tigre::worthPageLocking(isHostRegisterSupported, gpuData, splits, deviceCount)){  // as Siddon
        img_registered=tigre::hostRegister(img, (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*(size_t)sizeof(float));
        // The caller's array: left registered, the next call with it would fail.
        if (img_registered) onError.add([img]{ cudaHostUnregister(img); });
    }
    cudaCheckErrors("Error pinning memory");
#endif
    Point3D source, deltaU, deltaV, uvOrigin;
    // The geometry constants of each launch are filled here and copied to __constant__ memory before the kernel.
    // Ordinary host memory: page-locking these few KB (allocating and freeing them with cudaMallocHost and
    // cudaFreeHost) cost ~0.45 ms per call and did not change the kernels' time (measured on an RTX 4070: the
    // kernels read the __constant__ copy either way).
    
    Point3D* projParamsArrayHost = 0;
    projParamsArrayHost=(Point3D*)malloc(4*projPerLaunch*sizeof(Point3D));
    if (!projParamsArrayHost) mexErrMsgIdAndTxt("TIGRE:hostMemory","Not enough host memory for the geometry constants");
    onError.add([projParamsArrayHost]{ free(projParamsArrayHost); });
    float* projFloatsArrayHost = 0;
    projFloatsArrayHost=(float*)malloc(2*projPerLaunch*sizeof(float));
    if (!projFloatsArrayHost) mexErrMsgIdAndTxt("TIGRE:hostMemory","Not enough host memory for the geometry constants");
    onError.add([projFloatsArrayHost]{ free(projFloatsArrayHost); });
    cudaCheckErrors("Error allocating auxiliary constant memory");
    
    // Create Streams for overlapping memcopy and compute
    int nStream_device=2;
    int nStreams=deviceCount*nStream_device;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    memset(stream, 0, nStreams*sizeof(cudaStream_t));
    onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        for (int i = 0; i < nStream_device; ++i){
            cudaStreamCreate(&stream[i+dev*nStream_device]);
            
        }
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
        // They are all the same size, except the last one.
        linear_idx_start= (size_t)sp*(size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        CreateTextureInterp(gpuids,&img[linear_idx_start],geoArray[sp],d_cuArrTex,texImg,false,(gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice),stream);  // arrays allocated above
        cudaCheckErrors("Texture object creation fail");
        
        
        int divU,divV;
        tigre::detectorBlockShape(geo.nDetecU, geo.nDetecV, PIXEL_SIZE_BLOCK, &divU, &divV);  // thin detectors (fan beam)
        dim3 grid((geoArray[sp].nDetecU+divU-1)/divU,(geoArray[0].nDetecV+divV-1)/divV,projPerLaunch/PROJ_PER_BLOCK);
        dim3 block(divU,divV,PROJ_PER_BLOCK);  // PROJ_PER_BLOCK angles per block, projPerLaunch/PROJ_PER_BLOCK blocks in z
        
        unsigned int proj_global;
        float maxdist;
        // Now that we have prepared the image (piece of image) and parameters for kernels
        // we project for all angles.
        for (unsigned int i=0; i<noOfKernelCalls; i++) {
            for (dev=0;dev<deviceCount;dev++){
                float is_spherical=0;
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
                    
                    is_spherical+=abs(geoArray[sp].theta)+abs(geoArray[sp].psi);
                    
                    //precomute distances for faster execution
                    maxdist=maxdistanceCuboid(geoArray[sp],proj_global);
                    //Precompute per angle constant stuff for speed
                    computeDeltas(geoArray[sp], proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                    //Ray tracing!
                    projParamsArrayHost[4*j]=uvOrigin;		// 6*j because we have 6 Point3D values per projection
                    projParamsArrayHost[4*j+1]=deltaU;
                    projParamsArrayHost[4*j+2]=deltaV;
                    projParamsArrayHost[4*j+3]=source;
                    
                    projFloatsArrayHost[2*j]=geo.DSO[proj_global];
                    projFloatsArrayHost[2*j+1]=floor(maxdist);
                }
                
                cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*projPerLaunch,0,cudaMemcpyHostToDevice,stream[dev*nStream_device]);
                cudaMemcpyToSymbolAsync(projFloatsArrayDev, projFloatsArrayHost, sizeof(float)*2*projPerLaunch,0,cudaMemcpyHostToDevice,stream[dev*nStream_device]);
                cudaStreamSynchronize(stream[dev*nStream_device]);
                
                
                //TODO: we could do this around X and Y axis too, but we would need to compute the new axis of rotation (not possible to know from jsut the angles)
                if (!is_spherical){
                    kernelPixelDetector<false><<<grid,block,0,stream[dev*nStream_device]>>>(geoArray[sp],(directOut ? result[i*projPerLaunch] : dProjection[(i%2)+dev*2]),i,nangles_device,projPerLaunch,texImg[dev]);
                }
                else{
                    kernelPixelDetector<true> <<<grid,block,0,stream[dev*nStream_device]>>>(geoArray[sp],(directOut ? result[i*projPerLaunch] : dProjection[(i%2)+dev*2]),i,nangles_device,projPerLaunch,texImg[dev]);
                }
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
                    vecAddInPlaceInterp<<<(unsigned int)(((size_t)geo.nDetecU*geo.nDetecV*projection_this_block+MAXTREADS-1)/MAXTREADS),MAXTREADS,0,stream[dev*2]>>>(dProjection[(i%2)+dev*2],dProjection_accum[(i%2)+dev*2],(size_t)geo.nDetecU*geo.nDetecV*projection_this_block);
                }
            } // end accumulation case, where the image needs to be split 

            // Now, lets get out the projections from the previous execution of the kernels.
            if (i>0)
            {
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
                    if (!directOut) cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(i%2))+dev*2],  (size_t)projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost),stream[dev*2+copyStream]);
                }
            }
            // Make sure Computation on kernels has finished before we launch the next batch.
            // ...and the copy-out of the batch before it (copy stream), which reads the
            // buffer the next kernel writes. Today that copy goes to pageable host memory
            // and is already finished when cudaMemcpyAsync returns, but nothing else
            // guarantees it (a page-locked result would race), see #860.
            for (dev = 0; dev < deviceCount; dev++)
            {
                cudaSetDevice(gpuids[dev]);
                cudaStreamSynchronize(stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2+1]);
            }
        } // End noOfKernelCalls (i) loop.
        
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
            if (!directOut) cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(noOfKernelCalls%2))+dev*2], (size_t)projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost),stream[dev*2+copyStream]);
        }
        // Make sure everyone has done their bussiness before the next image split:
        for (dev = 0; dev < deviceCount; dev++)
        {
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
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
    }
    free(dProjection_accum);
    freeGeoArray(splits,geoArray);
    free(projParamsArrayHost);
    free(projFloatsArrayHost);
    
    
    for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]) ;
#ifndef NO_PINNED_MEMORY
    if (img_registered){
        cudaHostUnregister(img);
    }
#endif
    cudaCheckErrors("cudaFree  fail");
    
//     cudaDeviceReset();
    return 0;
}
// The copy goes on the stream of each GPU's kernels (stream[2*dev]), which orders it before them: the host does
// not wait for it.
void CreateTextureInterp(const GpuIds& gpuids,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,bool allocate, cudaMemcpyKind kind, cudaStream_t* stream)
{
    const unsigned int num_devices = gpuids.GetLength();
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    if(allocate){
        
        for (unsigned int dev = 0; dev < num_devices; dev++){
            cudaSetDevice(gpuids[dev]);
            
            //cudaArray Descriptor
            
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            //cuda Array
            tigre::checkMalloc(cudaMalloc3DArray(&d_cuArrTex[dev], &channelDesc, extent), "Ax:Interpolated_projection:GPUmemory", "Forward projection (interpolated)", (size_t)extent.width*extent.height*extent.depth*sizeof(float));
            cudaCheckErrors("Texture memory allocation fail");
        }
        
    }
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaMemcpy3DParms copyParams = {0};
        cudaSetDevice(gpuids[dev]);
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)imagedata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[dev];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
        cudaMemcpy3DAsync(&copyParams,stream[2*dev]);
        //cudaCheckErrors("Texture memory data copy fail");
        //Array creation End
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
        if (geo.accuracy>1){
            texDescr.filterMode = cudaFilterModePoint;
            geo.accuracy=1;
        }
        else{
            texDescr.filterMode = cudaFilterModeLinear;
        }
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
        cudaCheckErrors("Texture object creation fail");
    }
}

/* This code generates the geometries needed to split the image properly in
 * cases where the entire image does not fit in the memory of the GPU
 **/
void splitImageInterp(unsigned int splits,Geometry geo,Geometry* geoArray, unsigned int nangles){
    
    unsigned long splitsize=(geo.nVoxelZ+splits-1)/splits;// ceil if not divisible
    for(unsigned int sp=0;sp<splits;sp++){
        geoArray[sp]=geo;
        // All of them are splitsize, but the last one, possibly fewer (never none: see the callers). Signed: in
        // unsigned, nVoxelZ-splitsize*sp wraps around for a piece past the end.
        const long long rest=(long long)geo.nVoxelZ-(long long)splitsize*sp;
        geoArray[sp].nVoxelZ=(int)(rest<=0 ? 0 : (rest<(long long)splitsize ? rest : (long long)splitsize));
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
void computeDeltas(Geometry geo,unsigned int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){
    Point3D S;
    S.x=geo.DSO[i];
    S.y=0;
    S.z=0;
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD[i]-geo.DSO[i]);   P.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       P.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pu0.x=-(geo.DSD[i]-geo.DSO[i]);   Pu0.y= geo.dDetecU*(1-((float)geo.nDetecU/2)+0.5);       Pu0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pv0.x=-(geo.DSD[i]-geo.DSO[i]);   Pv0.y= geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       Pv0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-1);
    // Geomtric trasnformations:
    
    
    // Now we have the Real world (OXYZ) coordinates of the bottom corner and its two neighbours.
    // The obkjective is to get a position of the detector in a coordinate system where:
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
    P.x=0;Pu0.x=0;Pv0.x=0;
    
    // Roll pitch yaw
    rollPitchYaw(geo,i,&P);
    rollPitchYaw(geo,i,&Pu0);
    rollPitchYaw(geo,i,&Pv0);
    //Now ltes translate the detector coordinates to DOD (original position on real coordinate system:
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);
    Pu0.x=Pu0.x-(geo.DSD[i]-geo.DSO[i]);
    Pv0.x=Pv0.x-(geo.DSD[i]-geo.DSO[i]);
    //2: Offset detector
    
    
    //S doesnt need to chagne
    
    
    
    //3: Rotate around RZ RY RZ
    Point3D Pfinal, Pfinalu0, Pfinalv0;
    Pfinal.x  =P.x;
    Pfinal.y  =P.y  +geo.offDetecU[i]; Pfinal.z  =P.z  +geo.offDetecV[i];
    Pfinalu0.x=Pu0.x;
    Pfinalu0.y=Pu0.y  +geo.offDetecU[i]; Pfinalu0.z  =Pu0.z  +geo.offDetecV[i];
    Pfinalv0.x=Pv0.x;
    Pfinalv0.y=Pv0.y  +geo.offDetecU[i]; Pfinalv0.z  =Pv0.z  +geo.offDetecV[i];
    
    eulerZYZ(geo,&Pfinal);
    eulerZYZ(geo,&Pfinalu0);
    eulerZYZ(geo,&Pfinalv0);
    eulerZYZ(geo,&S);
    
    
    //3: Offset image (instead of offseting image, -offset everything else)
    
    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    Pfinalu0.x=Pfinalu0.x-geo.offOrigX[i];   Pfinalu0.y=Pfinalu0.y-geo.offOrigY[i];   Pfinalu0.z=Pfinalu0.z-geo.offOrigZ[i];
    Pfinalv0.x=Pfinalv0.x-geo.offOrigX[i];   Pfinalv0.y=Pfinalv0.y-geo.offOrigY[i];   Pfinalv0.z=Pfinalv0.z-geo.offOrigZ[i];
    S.x=S.x-geo.offOrigX[i];                 S.y=S.y-geo.offOrigY[i];                 S.z=S.z-geo.offOrigZ[i];
    
    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2-geo.dVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2-geo.dVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2-geo.dVoxelZ/2;
    Pfinalu0.x=Pfinalu0.x+geo.sVoxelX/2-geo.dVoxelX/2;    Pfinalu0.y=Pfinalu0.y+geo.sVoxelY/2-geo.dVoxelY/2;        Pfinalu0.z=Pfinalu0.z+geo.sVoxelZ/2-geo.dVoxelZ/2;
    Pfinalv0.x=Pfinalv0.x+geo.sVoxelX/2-geo.dVoxelX/2;    Pfinalv0.y=Pfinalv0.y+geo.sVoxelY/2-geo.dVoxelY/2;        Pfinalv0.z=Pfinalv0.z+geo.sVoxelZ/2-geo.dVoxelZ/2;
    S.x       =S.x+geo.sVoxelX/2-geo.dVoxelX/2;           S.y       =S.y+geo.sVoxelY/2-geo.dVoxelY/2;               S.z       =S.z      +geo.sVoxelZ/2-geo.dVoxelZ/2;
    
    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S.x       =S.x/geo.dVoxelX;           S.y       =S.y/geo.dVoxelY;             S.z       =S.z/geo.dVoxelZ;
    
    
    //mexPrintf("COR: %f \n",geo.COR[i]);
    //5. apply COR. Wherever everything was, now its offesetd by a bit.
//     Only wors for standard rotaiton, not aribtary axis rotation.
    float CORx, CORy;
    CORx=-geo.COR[i]*sin(geo.alpha)/geo.dVoxelX;
    CORy= geo.COR[i]*cos(geo.alpha)/geo.dVoxelY;
    Pfinal.x+=CORx;   Pfinal.y+=CORy;
    Pfinalu0.x+=CORx;   Pfinalu0.y+=CORy;
    Pfinalv0.x+=CORx;   Pfinalv0.y+=CORy;
    S.x+=CORx; S.y+=CORy;
    
    // return
    
    *uvorigin=Pfinal;
    
    deltaU->x=Pfinalu0.x-Pfinal.x;
    deltaU->y=Pfinalu0.y-Pfinal.y;
    deltaU->z=Pfinalu0.z-Pfinal.z;
    
    deltaV->x=Pfinalv0.x-Pfinal.x;
    deltaV->y=Pfinalv0.y-Pfinal.y;
    deltaV->z=Pfinalv0.z-Pfinal.z;
    
    *source=S;
}

float maxdistanceCuboid(Geometry geo,unsigned int i){
    ///////////
    // Compute initial "t" so we access safely as less as out of bounds as possible.
    //////////
    
    
    float maxCubX,maxCubY,maxCubZ;
    // Forgetting Z, compute mas distance: diagonal+offset
    maxCubX=(geo.nVoxelX/2+ abs(geo.offOrigX[i])/geo.dVoxelX);
    maxCubY=(geo.nVoxelY/2+ abs(geo.offOrigY[i])/geo.dVoxelY);
    maxCubZ=(geo.nVoxelZ/2+ abs(geo.offOrigZ[i])/geo.dVoxelZ);
    
    float a,b;
    a=geo.DSO[i]/geo.dVoxelX;
    b=geo.DSO[i]/geo.dVoxelY;
    
//  As the return of this value is in "voxel space", the source may have an elliptical curve.
//  The distance returned is the safe distance that can be skipped for a given angle alpha, before we need to start sampling.
    
    if (geo.theta==0.0f & geo.psi==0.0f) // Special case, it will make the code faster
        return max(a*b/sqrt(a*a*sin(geo.alpha)*sin(geo.alpha)+b*b*cos(geo.alpha)*cos(geo.alpha))-
                sqrt(maxCubX*maxCubX+maxCubY*maxCubY),0.0f);
    //TODO: think of more special cases?
    return max(geo.DSO[i]/max(max(geo.dVoxelX,geo.dVoxelY),geo.dVoxelZ)-sqrt(maxCubX*maxCubX+maxCubY*maxCubY+maxCubZ*maxCubZ),0.0f);

}
void rollPitchYaw(Geometry geo,unsigned int i, Point3D* point){
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
void eulerZYZ(Geometry geo,  Point3D* point){
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
    int deviceCount = gpuids.GetLength();
    for (int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) *mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            mexWarnMsgIdAndTxt("ray_interpolated_projection:ax:GPU","One (or more) of your GPUs is being heavily used by another program (less than half of its memory is free). TIGRE will try to run anyway, but it may need to split the problem into smaller pieces (slower) or run out of GPU memory.\n");
        }
        cudaCheckErrors("Check mem error");
        *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
    }
    *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
    
    //*mem_GPU_global= insert your known number here, in bytes.
}
