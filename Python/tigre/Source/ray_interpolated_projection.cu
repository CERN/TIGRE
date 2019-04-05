/*-------------------------------------------------------------------------
 *
 * CUDA functions for texture-memory interpolation based projection
 *
 * This file has the necesary fucntiosn to perform X-ray CBCT projection
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
#include "errors.hpp"
#include <math.h>
#include <stdio.h>

inline int cudaCheckErrors(const char * msg)
{
   cudaError_t __err = cudaGetLastError();
   if (__err != cudaSuccess)
   {
      printf("ray_interpolated_projection:%s:%s\n",msg, cudaGetErrorString(__err));
      cudaDeviceReset();
      return 1;
   }
   return 0;
}
    
    
// Declare the texture reference.
    
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
    void CreateTextureInterp(int num_devices,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,bool allocate);
__constant__ Point3D projParamsArrayDev[4*PROJ_PER_BLOCK];  // Dev means it is on device
__constant__ float projFloatsArrayDev[2*PROJ_PER_BLOCK];  // Dev means it is on device


__global__ void vecAddInPlaceInterp(float *a, float *b, unsigned long  n)
{
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    // Make sure we do not go out of bounds
    if (idx < n)
        a[idx] = a[idx] + b[idx];
}


template<bool sphericalrotation>
        __global__ void kernelPixelDetector( Geometry geo,
        float* detector,
        const int currProjSetNumber,
        const int totalNoOfProjections,
        cudaTextureObject_t tex){
    
    unsigned long  y = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long  x = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long projNumber=threadIdx.z;
    
    
    if ((x>= geo.nDetecU) | (y>= geo.nDetecV)|  (projNumber>=PROJ_PER_BLOCK))
        return;
    
    size_t idx =  (size_t)(x  * geo.nDetecV + y)+ (size_t)projNumber*geo.nDetecV *geo.nDetecU ;
    int indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    
    float DSO = projFloatsArrayDev[2*projNumber+0];
    float cropdist_init = projFloatsArrayDev[2*projNumber+1];
    
    
    
    /////// Get coordinates XYZ of pixel UV
    int pixelV = geo.nDetecV-y-1;
    int pixelU = x;
    
    
    
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
int interpolation_projection(float  *  img, Geometry geo, float** result,float const * const angles,int nangles){
    
    
    // Prepare for MultiGPU
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if(cudaCheckErrors("Device query fail")){return 1;}
    if (deviceCount == 0) {
        return ERR_NO_CAPABLE_DEVICES;
        }
    //
    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning trhown)
    int dev;
    char * devicenames;
    cudaDeviceProp deviceProp;
    
    for (dev = 0; dev < deviceCount; dev++) {
        cudaSetDevice(dev);
        cudaGetDeviceProperties(&deviceProp, dev);
        if (dev>0){
            if (strcmp(devicenames,deviceProp.name)!=0){
                printf("Ax:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed. \n Siddon_projection.cu line 275.");
                break;
            }
        }
        devicenames=deviceProp.name;
    }
    
    // Check free memory
    size_t mem_GPU_global;
    checkFreeMemory(deviceCount,&mem_GPU_global);
    
    size_t mem_image=(unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
    size_t mem_proj =(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV * sizeof(float);
    
    // Does everything fit in the GPUs?
    bool fits_in_memory=false;
    unsigned int splits=1;
    Geometry * geoArray;
    
    
    if (mem_image+2*PROJ_PER_BLOCK*mem_proj<mem_GPU_global){// yes it does
        fits_in_memory=true;
        geoArray=(Geometry*)malloc(sizeof(Geometry));
        geoArray[0]=geo;
    }
    else{// Nope nope.
        fits_in_memory=false; // Oh dear.
        // approx free memory we have. We already have left some extra 10% free for internal stuff
        // we need a second projection memory to combine multi-GPU stuff.
        size_t mem_free=mem_GPU_global-4*PROJ_PER_BLOCK*mem_proj;
        
        
        splits=mem_image/mem_free+1;// Ceil of the truncation
        geoArray=(Geometry*)malloc(splits*sizeof(Geometry));
        splitImageInterp(splits,geo,geoArray,nangles);
    }
    
    // Allocate auiliary memory for projections on the GPU to accumulate partial resutsl
    float ** dProjection_accum;
    size_t num_bytes_proj = PROJ_PER_BLOCK*geo.nDetecU*geo.nDetecV * sizeof(float);
    if (!fits_in_memory){
        dProjection_accum=(float**)malloc(2*deviceCount*sizeof(float*));
        for (dev = 0; dev < deviceCount; dev++) {
            cudaSetDevice(dev);
            for (int i = 0; i < 2; ++i){
                cudaMalloc((void**)&dProjection_accum[dev*2+i], num_bytes_proj);
                cudaMemset(dProjection_accum[dev*2+i],0,num_bytes_proj);
                if(cudaCheckErrors("cudaMallocauxiliarty projections fail")){return 1;}
            }
        }
    }
    
    // This is happening regarthless if the image fits on memory
    float** dProjection=(float**)malloc(2*deviceCount*sizeof(float*));
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        
        for (int i = 0; i < 2; ++i){
            cudaMalloc((void**)&dProjection[dev*2+i],   num_bytes_proj);
            cudaMemset(dProjection[dev*2+i]  ,0,num_bytes_proj);
            if(cudaCheckErrors("cudaMalloc projections fail")){return 1;}
        }
    }
    
    
    
    
    //Pagelock memory for syncronous copy.
    // Lets try to make the host memory pinned:
    // We laredy queried the GPU and assuemd they are the same, thus shoudl have the same attributes.
    int isHostRegisterSupported;
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,0);
    // empirical testing shows that when the image split is smaller than 1 (also implies the image is not very big), the time to
    // pin the memory is greater than the lost time in Syncronously launching the memcpys. This is only worth it when the image is too big.
    if (isHostRegisterSupported & splits>1){
        cudaHostRegister(img, (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*(size_t)sizeof(float),cudaHostRegisterPortable);
    }
    
    
    
    
    Point3D source, deltaU, deltaV, uvOrigin;
    
    Point3D* projParamsArrayHost;
    cudaMallocHost((void**)&projParamsArrayHost,4*PROJ_PER_BLOCK*sizeof(Point3D));
    float* projFloatsArrayHost;
    cudaMallocHost((void**)&projFloatsArrayHost,2*PROJ_PER_BLOCK*sizeof(float));
    
    
    // Create Streams for overlapping memcopy and compute
    int nStream_device=2;
    int nStreams=deviceCount*nStream_device;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        for (int i = 0; i < nStream_device; ++i){
            cudaStreamCreate(&stream[i+dev*nStream_device]);
            
        }
    }
    if(cudaCheckErrors("Stream creation fail")){return 1;}
    int nangles_device=(nangles+deviceCount-1)/deviceCount;
    int nangles_last_device=(nangles-(deviceCount-1)*nangles_device);
    unsigned int noOfKernelCalls = (nangles_device+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_BLOCK
    unsigned int last_device_blocks= (nangles_last_device+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK; // we will use this in the memory management.
    int projection_this_block;

    
    cudaTextureObject_t *texImg = new cudaTextureObject_t[deviceCount];
    cudaArray **d_cuArrTex = new cudaArray*[deviceCount];
    for (unsigned int sp=0;sp<splits;sp++){
        
        // Create texture objects for all GPUs
        
        
        size_t linear_idx_start;
        //First one shoudl always be  the same size as all the rest but the last
        linear_idx_start= (size_t)sp*(size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        CreateTextureInterp(deviceCount,&img[linear_idx_start],geoArray[sp],d_cuArrTex,texImg,!sp);
        if(cudaCheckErrors("Texture object creation fail")){return 1;}
        
        
        int divU,divV;
        divU=PIXEL_SIZE_BLOCK;
        divV=PIXEL_SIZE_BLOCK;
        dim3 grid((geoArray[sp].nDetecU+divU-1)/divU,(geoArray[0].nDetecV+divV-1)/divV,1);
        dim3 block(divU,divV,PROJ_PER_BLOCK);
        
        unsigned int proj_global;
        unsigned int i;
        float maxdist;
        // Now that we have prepared the image (piece of image) and parameters for kernels
        // we project for all angles.
        for ( i=0; i<noOfKernelCalls; i++){
            for (dev=0;dev<deviceCount;dev++){
                float is_spherical=0;
                cudaSetDevice(dev);
                
                for(unsigned int j=0; j<PROJ_PER_BLOCK; j++){
                    proj_global=(i*PROJ_PER_BLOCK+j)+dev*nangles_device;
                    if (proj_global>=nangles)
                        break;
                     if ((i*PROJ_PER_BLOCK+j)>=nangles_device)
                        break;
                    geo.alpha=angles[proj_global*3];
                    geo.theta=angles[proj_global*3+1];
                    geo.psi  =angles[proj_global*3+2];
                    
                    is_spherical+=abs(geo.theta)+abs(geo.psi);
                    
                    //precomute distances for faster execution
                    maxdist=maxdistanceCuboid(geo,proj_global);
                    //Precompute per angle constant stuff for speed
                    computeDeltas(geo, proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                    //Ray tracing!
                    projParamsArrayHost[4*j]=uvOrigin;		// 6*j because we have 6 Point3D values per projection
                    projParamsArrayHost[4*j+1]=deltaU;
                    projParamsArrayHost[4*j+2]=deltaV;
                    projParamsArrayHost[4*j+3]=source;
                    
                    projFloatsArrayHost[2*j]=geo.DSO[proj_global];
                    projFloatsArrayHost[2*j+1]=floor(maxdist);
                }
                
                cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[dev*nStream_device]);
                cudaMemcpyToSymbolAsync(projFloatsArrayDev, projFloatsArrayHost, sizeof(float)*2*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[dev*nStream_device]);
                cudaStreamSynchronize(stream[dev*nStream_device]);
                
                
                //TODO: we could do this around X and Y axis too, but we would need to compute the new axis of rotation (not possible to know from jsut the angles)
                if (!is_spherical){
                    kernelPixelDetector<false><<<grid,block,0,stream[dev*nStream_device]>>>(geoArray[sp],dProjection[(i%2)+dev*2],i,nangles_device,texImg[dev]);
                }
                else{
                    kernelPixelDetector<true> <<<grid,block,0,stream[dev*nStream_device]>>>(geoArray[sp],dProjection[(i%2)+dev*2],i,nangles_device,texImg[dev]);
                }
            }
            
            
            // Now that the computation is happening, we need to either prepare the memory for
            // combining of the projections (splits>1) or start removing previous results.
            
            
            
            // If our image does not fit in memory then we need to make sure we accumulate previous results too.
            if( !fits_in_memory&&sp>0){
                // First, grab previous results and put them in the auxiliary variable
                for (dev = 0; dev < deviceCount; dev++){
                    projection_this_block=PROJ_PER_BLOCK;
                    cudaSetDevice(dev);
                    // this werid code makes sure we dont access bad memory. Its necesary for deviceCount>2
                    if (dev+1==deviceCount){ // if its the last device
                        if(i+1==last_device_blocks) // If we are in the last block of the last device, how many projections?
                            projection_this_block=nangles_last_device-(last_device_blocks-1)*PROJ_PER_BLOCK;
                        if(i+1>last_device_blocks) // As the last device can have less blocs, i may be over it.
                            break;
                    }else{
                        if(i+1==noOfKernelCalls) // if its not the last device, it can still be the lat block
                            projection_this_block=nangles_device-(noOfKernelCalls-1)*PROJ_PER_BLOCK;
                    }
                    cudaMemcpyAsync(dProjection_accum[(i%2)+dev*2], result[i*PROJ_PER_BLOCK+dev*nangles_device], projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), cudaMemcpyHostToDevice,stream[dev*2+1]);
                }
                // Second, take the results from current compute call and add it to the code in execution.
                for (dev = 0; dev < deviceCount; dev++){
                    
                    projection_this_block=PROJ_PER_BLOCK;
                    cudaSetDevice(dev);
                    // this werid code makes sure we dont access bad memory. Its necesary for deviceCount>2
                    if (dev+1==deviceCount){ // if its the last device
                        if(i+1==last_device_blocks) // If we are in the last block of the last device, how many projections?
                            projection_this_block=nangles_last_device-(last_device_blocks-1)*PROJ_PER_BLOCK;
                        if(i+1>last_device_blocks) // As the last device can have less blocs, i may be over it.
                            break;
                    }else{
                        if(i+1==noOfKernelCalls) // if its not the last device, it can still be the lat block
                            projection_this_block=nangles_device-(noOfKernelCalls-1)*PROJ_PER_BLOCK;
                    }
                    cudaStreamSynchronize(stream[dev*2+1]); // wait until copy is finished
                    vecAddInPlaceInterp<<<(geo.nDetecU*geo.nDetecV*projection_this_block+MAXTREADS-1)/MAXTREADS,MAXTREADS,0,stream[dev*2]>>>(dProjection[(i%2)+dev*2],dProjection_accum[(i%2)+dev*2],(unsigned long)geo.nDetecU*geo.nDetecV*projection_this_block);
                }
            }
            // Now, lets get out the projections from the previous execution of the kernels.
            if (i>0){
                for (dev = 0; dev < deviceCount; dev++){
                    projection_this_block=PROJ_PER_BLOCK;
                    cudaSetDevice(dev);
                    if (dev+1==deviceCount && i+1==noOfKernelCalls && last_device_blocks!=noOfKernelCalls){ 
                            projection_this_block=nangles_last_device-(last_device_blocks-1)*PROJ_PER_BLOCK;
                    }
                    cudaMemcpyAsync(result[(i-1)*PROJ_PER_BLOCK+dev*nangles_device], dProjection[(int)(!(i%2))+dev*2],  projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*2+1]);
                }
            }
            // Make sure Computation on kernels has finished before we launch the next batch.
            for (dev = 0; dev < deviceCount; dev++){
                cudaSetDevice(dev);
                cudaStreamSynchronize(stream[dev*2]);
            }
        }
        
        // We still have the last one to get out, do that one
        int projection_this_block;
        for (dev = 0; dev < deviceCount; dev++){
            projection_this_block=PROJ_PER_BLOCK;
            cudaSetDevice(dev);
            // this werid code makes sure we dont access bad memory. Its necesary for deviceCount>2
            if (dev+1==deviceCount){ // if its the last device
                projection_this_block=nangles_last_device-(last_device_blocks-1)*PROJ_PER_BLOCK;
                if(i>last_device_blocks) // As the last device can have less blocs, i may be over it.
                    break;
            }else{
               projection_this_block=nangles_device-(noOfKernelCalls-1)*PROJ_PER_BLOCK;
            }
            cudaDeviceSynchronize();
            if(cudaCheckErrors("Fail memcopy fail")){return 1;}
            cudaMemcpyAsync(result[(i-1)*PROJ_PER_BLOCK+dev*nangles_device], dProjection[(int)(!(i%2))+dev*2], projection_this_block*geo.nDetecV*geo.nDetecU*sizeof(float), cudaMemcpyDeviceToHost,stream[dev*2+1]);
        }
        // Free memory for the next piece of image
        
        cudaDeviceSynchronize();
        
    }
    
    if(cudaCheckErrors("Main loop  fail")){return 1;}
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaDestroyTextureObject(texImg[dev]);
        cudaFreeArray(d_cuArrTex[dev]);
    }
    // Freeing Stage
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaFree(dProjection[dev*2]);
        cudaFree(dProjection[dev*2+1]);
        
    }
    free(dProjection);
    
    if(!fits_in_memory){
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(dev);
            cudaFree(dProjection_accum[dev*2]);
            cudaFree(dProjection_accum[dev*2+1]);
            
        }
        free(dProjection_accum);
    }
    freeGeoArray(splits,geoArray);
    cudaFreeHost(projParamsArrayHost);
    
    
    for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]) ;
    
    if (isHostRegisterSupported & splits>1){
        cudaHostUnregister(img);
    }
    if(cudaCheckErrors("cudaFree  fail")){return 1;}
    
//     cudaDeviceReset();
    return 0;
}
void CreateTextureInterp(int num_devices,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,bool allocate)
{
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    if(allocate){
        
        for (unsigned int i = 0; i < num_devices; i++){
            cudaSetDevice(i);
            
            //cudaArray Descriptor
            
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            //cuda Array
            cudaMalloc3DArray(&d_cuArrTex[i], &channelDesc, extent);
            cudaCheckErrors("Texture memory allocation fail");
        }
        
    }
    for (unsigned int i = 0; i < num_devices; i++){
        cudaMemcpy3DParms copyParams = {0};
        cudaSetDevice(i);
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)imagedata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[i];
        copyParams.extent   = extent;
        copyParams.kind     = cudaMemcpyHostToDevice;
        cudaMemcpy3DAsync(&copyParams);
        //cudaCheckErrors("Texture memory data copy fail");
        //Array creation End
    }
    for (unsigned int i = 0; i < num_devices; i++){
        cudaSetDevice(i);
        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArrTex[i];
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
        cudaCreateTextureObject(&texImage[i], &texRes, &texDescr, NULL);
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
void checkFreeMemory(int deviceCount,size_t *mem_GPU_global){
    size_t memfree;
    size_t memtotal;
    
    for (int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) *mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            printf("tvDenoise:tvdenoising:GPU","One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n Free the GPU to run TIGRE\n");
        }
        cudaCheckErrors("Check mem error");
        *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
    }
    *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
    
    //*mem_GPU_global= insert your known number here, in bytes.
}