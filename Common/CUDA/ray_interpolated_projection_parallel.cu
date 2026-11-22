/*-------------------------------------------------------------------------
 *
 * CUDA functions for texture-memory interpolation based projection
 *
 * This file has the necessary functions to perform X-ray parallel projection
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
#include <vector>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "ray_interpolated_projection_parallel.hpp"
#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include "launch_shapes.hpp"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                tigreCudaError("TIGRE:Ax:interpolated_parallel", msg, cudaGetErrorString(__err));\
        } \
} while (0)
    
    

#define MAXTREADS 1024
#define PROJ_PER_BLOCK 8
#define PIXEL_SIZE_BLOCK 8
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
void CreateTextureParallelInterp(const int* devs, int ndev, float* image, Geometry geo, cudaArray** d_cuArrTex, cudaTextureObject_t *texImage, cudaStream_t* stream, cudaMemcpyKind kind);
__constant__ Point3D projParamsArrayDev[4*PROJ_PER_BLOCK];  // Dev means it is on device
__constant__ float projFloatsArrayDev[2*PROJ_PER_BLOCK];  // Dev means it is on device



__global__ void kernelPixelDetector_parallel_interpolated( Geometry geo,
        float* detector,
        const int currProjSetNumber, const int totalNoOfProjections, cudaTextureObject_t tex)
{
//         Point3D source ,
//         Point3D deltaU,
//         Point3D deltaV,
//         Point3D uvOrigin,
//         float DSO,
//         float maxdist){
    
    unsigned long long u = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long v = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long long projNumber=threadIdx.z;
    
    if (u>= geo.nDetecU || v>= geo.nDetecV || projNumber>=PROJ_PER_BLOCK)
        return;
    
    int indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    
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
    
    float DSO = projFloatsArrayDev[2*projNumber+0];
    float maxdist = projFloatsArrayDev[2*projNumber+1];
    
    
    /////// Get coordinates XYZ of pixel UV
    unsigned long pixelV = geo.nDetecV-v-1;
    unsigned long pixelU = u;
    
    
    float vectX,vectY,vectZ;
    Point3D P;
    P.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    P.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    P.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
    Point3D S;
    S.x=(source.x+pixelU*deltaU.x+pixelV*deltaV.x);
    S.y=(source.y+pixelU*deltaU.y+pixelV*deltaV.y);
    S.z=(source.z+pixelU*deltaU.z+pixelV*deltaV.z);
    
    // Length is the ray length in normalized space
    double length=sqrtf((S.x-P.x)*(S.x-P.x)+(S.y-P.y)*(S.y-P.y)+(S.z-P.z)*(S.z-P.z));
    //now legth is an integer of Nsamples that are required on this line
    length=ceilf(length/geo.accuracy);//Divide the directional vector by an integer
    vectX=(P.x -S.x)/(length);
    vectY=(P.y -S.y)/(length);
    vectZ=(P.z -S.z)/(length);
    
    
//     //Integrate over the line
    float tx,ty,tz;
    float sum=0;
    float i;
    
    
    // limit the amount of mem access after the cube, but before the detector.
    if ((2*DSO/geo.dVoxelX+maxdist)/geo.accuracy  <   length)
        length=ceilf((2*DSO/geo.dVoxelX+maxdist)/geo.accuracy);
    //Length is not actually a length, but the amount of memreads with given accuracy ("samples per voxel")
    
    for (i=floorf(maxdist/geo.accuracy); i<=length; i=i+1){
        tx=vectX*i+S.x;
        ty=vectY*i+S.y;
        tz=vectZ*i+S.z;
        
        sum += tex3D<float>(tex, tx+0.5f, ty+0.5f, tz+0.5f); // this line is 94% of time.
        
    }
    float deltalength=sqrtf((vectX*geo.dVoxelX)*(vectX*geo.dVoxelX)+
            (vectY*geo.dVoxelY)*(vectY*geo.dVoxelY)+
            (vectZ*geo.dVoxelZ)*(vectZ*geo.dVoxelZ) );
    detector[idx]=sum*deltalength;
}



// Adds b to a, element by element (the partial projections of the image pieces).
static __global__ void vecAddInPlaceParallelInterp(float *a, const float *b, size_t n)
{
    size_t idx = (size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if (idx < n)
        a[idx] = a[idx] + b[idx];
}

// Free memory on the most used of these GPUs, less 5% for CUDA's own needs.
static size_t freeMemoryParallelInterp(const int* devs, int ndev){
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

int interpolation_projection_parallel(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // Same scheme as the cone-beam projector (ray_interpolated_projection.cu): the angles are shared between the GPUs, each GPU
    // holds the whole image (or, if it does not fit, one z-piece of it at a time, adding up the projections of the
    // pieces), and two projection buffers per GPU let the copy-out of one batch overlap the next (#286).
    int deviceCount = gpuids.GetLength();
    // If a CUDA error raises (Python build) before the end of this function,
    // onError releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Interpolated_projection_parallel:GPUselect","There are no available device(s) that support CUDA\n");
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
        mexWarnMsgIdAndTxt("Ax:Interpolated_projection_parallel:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;

    // Does the image fit on each GPU, with the projection buffers? If not, split it in z.
    const size_t mem_GPU_global=freeMemoryParallelInterp(devIds.data(), deviceCount);
    const size_t mem_image=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*sizeof(float);
    const size_t num_bytes_proj=(size_t)PROJ_PER_BLOCK*(size_t)geo.nDetecU*(size_t)geo.nDetecV*sizeof(float);
    unsigned int splits=1;
    if (mem_image+2*num_bytes_proj>=mem_GPU_global){
        // Split: 2 more projection buffers to add up the pieces.
        const size_t mem_slice=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*sizeof(float);
        if (mem_GPU_global <= 4*num_bytes_proj+mem_slice)
            tigreNotEnoughGpuMemory("Ax:Interpolated_projection_parallel:GPUmemory","Forward projection (interpolated, parallel beam)",4*num_bytes_proj+mem_slice,mem_GPU_global);
        splits=mem_image/(mem_GPU_global-4*num_bytes_proj)+1;
        if (splits>geo.nVoxelZ) splits=geo.nVoxelZ;  // at least one slice per piece
    }
    // The image texture (a 3D array on each GPU, for the largest piece) is allocated first, with the projection
    // buffers tried on top. A 3D array needs more GPU memory than its raw size, so it can fail where the estimate
    // above said it fits: then split the image in more pieces and try again, down to one slice per piece.
    cudaArray **d_cuArrTex = new cudaArray*[deviceCount]();
    Geometry* geoArray = NULL;
    while (true){
        geoArray = (Geometry*)malloc(splits*sizeof(Geometry));
        splitImageInterp(splits,geo,geoArray,nangles);
        const cudaExtent extent = make_cudaExtent(geoArray[0].nVoxelX, geoArray[0].nVoxelY, geoArray[0].nVoxelZ);
        if (tigre::tryMalloc3DArrays(devIds.data(), deviceCount, d_cuArrTex, extent, (splits>1 ? 4 : 2), num_bytes_proj, "Ax:Interpolated_projection_parallel:GPU"))
            break;
        const unsigned long piece=geoArray[0].nVoxelZ;
        freeGeoArray(splits, geoArray);
        if (piece<=1 || splits>=geo.nVoxelZ){
            delete[] d_cuArrTex;
            size_t free_bytes=0, total_bytes=0;
            cudaMemGetInfo(&free_bytes, &total_bytes);
            tigreNotEnoughGpuMemory("Ax:Interpolated_projection_parallel:GPUmemory","Forward projection (interpolated, parallel beam)",4*num_bytes_proj+(size_t)geo.nVoxelX*geo.nVoxelY*sizeof(float),free_bytes);
        }
        while (splits<geo.nVoxelZ && (geo.nVoxelZ+splits-1)/splits>=piece) splits++;  // smaller pieces
    }
    const bool fits_in_memory = (splits==1);
    onError.add([geoArray, splits]{ freeGeoArray(splits, geoArray); });
    onError.add([d_cuArrTex, deviceCount, devIds]{
        for (int d = 0; d < deviceCount; d++) { cudaSetDevice(devIds[d]); if (d_cuArrTex[d]) cudaFreeArray(d_cuArrTex[d]); }
        delete[] d_cuArrTex;
    });

    // Buffers to add up the projections of the image pieces.
    float ** dProjection_accum = 0;
    if (!fits_in_memory){
        dProjection_accum=(float**)malloc(2*deviceCount*sizeof(float*));
        memset(dProjection_accum, 0, 2*deviceCount*sizeof(float*));
        onError.add([dProjection_accum, deviceCount, devIds]{
            for (int s = 0; s < 2*deviceCount; s++) { cudaSetDevice(devIds[s/2]); cudaFree(dProjection_accum[s]); }
            free(dProjection_accum);
        });
        for (dev = 0; dev < deviceCount; dev++) {
            cudaSetDevice(devIds[dev]);
            for (int i = 0; i < 2; ++i){
                tigre::checkMalloc(cudaMalloc((void**)&dProjection_accum[dev*2+i], num_bytes_proj), "Ax:Interpolated_projection_parallel:GPUmemory", "Forward projection (interpolated, parallel beam)", num_bytes_proj);
                cudaMemset(dProjection_accum[dev*2+i],0,num_bytes_proj); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
            }
        }
        cudaCheckErrors("cudaMalloc auxiliary projections fail");
    }

    float** dProjection=(float**)malloc(2*deviceCount*sizeof(float*));
    memset(dProjection, 0, 2*deviceCount*sizeof(float*));
    onError.add([dProjection, deviceCount, devIds]{
        for (int s = 0; s < 2*deviceCount; s++) { cudaSetDevice(devIds[s/2]); cudaFree(dProjection[s]); }
        free(dProjection);
    });
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(devIds[dev]);
        for (int i = 0; i < 2; ++i){
            tigre::checkMalloc(cudaMalloc((void**)&dProjection[dev*2+i], num_bytes_proj), "Ax:Interpolated_projection_parallel:GPUmemory", "Forward projection (interpolated, parallel beam)", num_bytes_proj);
            cudaMemset(dProjection[dev*2+i],0,num_bytes_proj); cudaDeviceSynchronize(); // built with --default-stream=per-thread: the memset is not ordered with our streams
        }
    }
    cudaCheckErrors("cudaMalloc projections fail");

    // Page-lock the image when it is copied in pieces or to several GPUs: the copies are then faster and asynchronous.
    bool registered=false;
#ifndef NO_PINNED_MEMORY
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,devIds[0]);
#endif
    if (!gpuData && isHostRegisterSupported && (splits>1 || deviceCount>1)){
        registered=tigre::hostRegister(img, mem_image);
        // The caller's array: left registered, the next call with it would fail.
        if (registered) onError.add([img]{ cudaHostUnregister(img); });
    }
#endif
    cudaCheckErrors("Error pinning memory");

    Point3D source, deltaU, deltaV, uvOrigin;
    Point3D* projParamsArrayHost = 0;
    cudaMallocHost((void**)&projParamsArrayHost,4*PROJ_PER_BLOCK*sizeof(Point3D));
    onError.add([projParamsArrayHost]{ cudaFreeHost(projParamsArrayHost); });
    float* projFloatsArrayHost = 0;
    cudaMallocHost((void**)&projFloatsArrayHost,2*PROJ_PER_BLOCK*sizeof(float));
    onError.add([projFloatsArrayHost]{ cudaFreeHost(projFloatsArrayHost); });
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

    // Value-initialised (zeroed) so CreateTextureParallelInterp can tell an empty slot from one
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
        CreateTextureParallelInterp(devIds.data(),deviceCount,&img[linear_idx_start],geoArray[sp],d_cuArrTex,texImg,stream,kindIn);
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
                    const float maxdist=maxdistanceCuboid(geoArray[sp],proj_global);
                    computeDeltas_parallel(geoArray[sp],geoArray[sp].alpha,proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                    projParamsArrayHost[4*j]=uvOrigin;
                    projParamsArrayHost[4*j+1]=deltaU;
                    projParamsArrayHost[4*j+2]=deltaV;
                    projParamsArrayHost[4*j+3]=source;
                    projFloatsArrayHost[2*j]=geoArray[sp].DSO[proj_global];
                    projFloatsArrayHost[2*j+1]=floor(maxdist);
                }
                cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[dev*2]);
                cudaMemcpyToSymbolAsync(projFloatsArrayDev, projFloatsArrayHost, sizeof(float)*2*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2]);  // the host arrays are refilled for the next GPU
                kernelPixelDetector_parallel_interpolated<<<numBlocks,threadsPerBlock,0,stream[dev*2]>>>(geoArray[sp],dProjection[(i%2)+dev*2],i,nangles_device,texImg[dev]);
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
                    vecAddInPlaceParallelInterp<<<(unsigned int)((len+MAXTREADS-1)/MAXTREADS),MAXTREADS,0,stream[dev*2]>>>(dProjection[(i%2)+dev*2],dProjection_accum[(i%2)+dev*2],len);
                }
            }

            // Copy out the previous batch while this one computes.
            if (i>0){
                for (dev = 0; dev < deviceCount; dev++){
                    const int n=projectionsInBatch(i-1,dev);
                    if (n<0) break;
                    cudaSetDevice(devIds[dev]);
                    const unsigned int proj_global=(i-1)*PROJ_PER_BLOCK+dev*nangles_device;
                    cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(i%2))+dev*2], (size_t)n*geo.nDetecV*geo.nDetecU*sizeof(float), kindOut, stream[dev*2+copyStream]);
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
            cudaMemcpyAsync(result[proj_global], dProjection[(int)(!(noOfKernelCalls%2))+dev*2], (size_t)n*geo.nDetecV*geo.nDetecU*sizeof(float), kindOut, stream[dev*2+copyStream]);
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
    if (!fits_in_memory) free(dProjection_accum);
    free(stream);
    freeGeoArray(splits,geoArray);
    cudaFreeHost(projParamsArrayHost);
    cudaFreeHost(projFloatsArrayHost);
    if (registered) cudaHostUnregister(img);
    cudaCheckErrors("cudaFree fail");
    return 0;
}




/* This code precomputes The location of the source and the Delta U and delta V (in the warped space)
 * to compute the locations of the x-rays. While it seems verbose and overly-optimized,
 * it does saves about 30% of each of the kernel calls. Thats something!
 **/
void computeDeltas_parallel(Geometry geo, float alpha,unsigned int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){
    Point3D S;
    S.x=geo.DSO[i];
    S.y=geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);
    S.z=geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    
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
    
    
    //3: Rotate (around z)!
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
    
    
    
    //2: Offset image (instead of offseting image, -offset everything else)
    
    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    Pfinalu0.x=Pfinalu0.x-geo.offOrigX[i];   Pfinalu0.y=Pfinalu0.y-geo.offOrigY[i];   Pfinalu0.z=Pfinalu0.z-geo.offOrigZ[i];
    Pfinalv0.x=Pfinalv0.x-geo.offOrigX[i];   Pfinalv0.y=Pfinalv0.y-geo.offOrigY[i];   Pfinalv0.z=Pfinalv0.z-geo.offOrigZ[i];
    S.x=S.x-geo.offOrigX[i];       S.y=S.y-geo.offOrigY[i];       S.z=S.z-geo.offOrigZ[i];
    
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
    
    
    
    //5. apply COR. Wherever everything was, now its offesetd by a bit
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
void CreateTextureParallelInterp(const int* devs, int ndev, float* image, Geometry geo, cudaArray** d_cuArrTex, cudaTextureObject_t *texImage, cudaStream_t* stream, cudaMemcpyKind kind){
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
            tigre::checkMalloc(cudaMalloc3DArray(&d_cuArrTex[dev], &channelDesc, extent), "Ax:Interpolated_projection_parallel:GPUmemory", "Forward projection (interpolated, parallel beam)", (size_t)extent.width*extent.height*extent.depth*sizeof(float));
        }
        cudaMemcpy3DParms copyParams = {0};
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)image, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[dev];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
        cudaMemcpy3DAsync(&copyParams,stream[dev*2+1]);
    }
    for (int dev = 0; dev < ndev; dev++){
        cudaSetDevice(devs[dev]);
        // Wait for the copy: the kernels run on the other stream and would otherwise read the texture while it is
        // still being filled (always for GPU arrays, whose device-to-device copies are fully asynchronous, and for
        // pinned host memory).
        cudaStreamSynchronize(stream[dev*2+1]);
        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArrTex[dev];
        cudaTextureDesc     texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = false;
        texDescr.filterMode = cudaFilterModeLinear;
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
