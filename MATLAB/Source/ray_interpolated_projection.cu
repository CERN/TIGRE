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
#include "mex.h"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                cudaDeviceReset();\
                mexErrMsgIdAndTxt("TIGRE:Ax:interpolated",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
// Declare the texture reference.

#define MAXTREADS 1024
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
void CreateTextureInterp(int num_devices,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage);


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
        Point3D source ,
        Point3D deltaU,
        Point3D deltaV,
        Point3D uvOrigin,
        float DSO,
        float cropdist_init,
        cudaTextureObject_t tex){
    
    unsigned long  y = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long  x = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long  idx =  x  * geo.nDetecV + y;
    
    if ((x>= geo.nDetecU) | (y>= geo.nDetecV))
        return;
    
    
    
    
    /////// Get coordinates XYZ of pixel UV
    int pixelV = geo.nDetecV-y-1;
    int pixelU = x;
    
    
    
    float vectX,vectY,vectZ;
    Point3D P;
    P.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    P.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    P.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
    
    // Length is the ray length in normalized space
    float length=sqrtf((source.x-P.x)*(source.x-P.x)+(source.y-P.y)*(source.y-P.y)+(source.z-P.z)*(source.z-P.z));
    //now legth is an integer of Nsamples that are required on this line
    length=ceilf(length/geo.accuracy);//Divide the directional vector by an integer
    vectX=(P.x -source.x)/(length);
    vectY=(P.y -source.y)/(length);
    vectZ=(P.z -source.z)/(length);
    
    
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
int interpolation_projection(float const * const img, Geometry geo, float** result,float const * const angles,int nangles){
    
    
    // Prepare for MultiGPU
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    cudaCheckErrors("Device query fail");
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Interpolated_projection:GPUselect","There are no available device(s) that support CUDA\n");
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
                mexWarnMsgIdAndTxt("Ax:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed. \n Siddon_projection.cu line 275.");
                break;
            }
        }
        devicenames=deviceProp.name;
    }
    
    // Check free memory
    size_t memfree;
    size_t memtotal;
    size_t mem_GPU_global;
    for (unsigned int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            mexErrMsgIdAndTxt("Ax:GPUselect","One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n Free the GPU to run TIGRE\n");
        }
        cudaCheckErrors("Check mem error");
        
        mem_GPU_global=(memfree<mem_GPU_global)?memfree:mem_GPU_global;
    }
    mem_GPU_global=(size_t)((double)mem_GPU_global*0.95); // lets leave 10% for the GPU. Too much? maybe, but probably worth saving.

    size_t mem_image=(unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
    size_t mem_proj =(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV * sizeof(float);
    
    // Does everything fit in the GPUs?
    bool fits_in_memory=false;
    unsigned int splits=1;
    Geometry * geoArray;
    
  
    if (mem_image+2*mem_proj<mem_GPU_global){// yes it does
        fits_in_memory=true;
        geoArray=(Geometry*)malloc(sizeof(Geometry));
        geoArray[0]=geo;
    }
    else{// Nope nope.
        fits_in_memory=false; // Oh dear.
        // approx free memory we have. We already have left some extra 10% free for internal stuff
        // we need a second projection memory to combine multi-GPU stuff.
        size_t mem_free=mem_GPU_global-4*mem_proj;
        

        splits=mem_image/mem_free+1;// Ceil of the truncation
        geoArray=(Geometry*)malloc(splits*sizeof(Geometry));
        splitImageInterp(splits,geo,geoArray,nangles);
    }
    
    // Allocate auiliary memory for projections on the GPU to accumulate partial resutsl
    float ** dProjection_accum;
    size_t num_bytes_proj = geo.nDetecU*geo.nDetecV * sizeof(float);
    if (!fits_in_memory){
        dProjection_accum=(float**)malloc(2*deviceCount*sizeof(float*));
        for (dev = 0; dev < deviceCount; dev++) {
            cudaSetDevice(dev);
            for (int i = 0; i < 2; ++i){
                cudaMalloc((void**)&dProjection_accum[dev*deviceCount+i], num_bytes_proj);
                cudaMemset(dProjection_accum[dev*deviceCount+i],0,num_bytes_proj);
                cudaCheckErrors("cudaMallocauxiliarty projections fail");
            }
        }
    }
    
    // This is happening regarthless if the image fits on memory
    float** dProjection=(float**)malloc(2*deviceCount*sizeof(float*));
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        
        for (int i = 0; i < 2; ++i){
            cudaMalloc((void**)&dProjection[dev*deviceCount+i],   num_bytes_proj);
            cudaMemset(dProjection[dev*deviceCount+i]  ,0,num_bytes_proj);
            cudaCheckErrors("cudaMalloc projections fail");
        }
    }
    
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
    cudaCheckErrors("Stream creation fail");
    
    
    float cropdist_init;
    for (unsigned int sp=0;sp<splits;sp++){
        
        // Create texture objects for all GPUs
        cudaTextureObject_t *texImg = new cudaTextureObject_t[deviceCount];
        cudaArray **d_cuArrTex = new cudaArray*[deviceCount];
        
        size_t linear_idx_start;
        //First one shoudl always be  the same size as all the rest but the last
        linear_idx_start= (size_t)sp*(size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        CreateTextureInterp(deviceCount,&img[linear_idx_start],geoArray[sp],d_cuArrTex,texImg);
        cudaCheckErrors("Texture object creation fail");
        
        
        
        
        Point3D source, deltaU, deltaV, uvOrigin;
        
        
        int divU,divV;
        divU=8;
        divV=8;
        dim3 grid((geoArray[sp].nDetecU+divU-1)/divU,(geoArray[0].nDetecV+divV-1)/divV,1);
        dim3 block(divU,divV,1);
        
        
        
        // Now that we have prepared the image (piece of image) and parameters for kernels
        // we project for all angles.
        for (unsigned int i=0;i<nangles;i+=(unsigned int)deviceCount){
            
            // First, we launch all the kernels for this angle

            for (dev = 0; dev < deviceCount; dev++){
                if(i+dev<nangles){
                    geoArray[sp].alpha=angles[(i+dev)*3];
                    geoArray[sp].theta=angles[(i+dev)*3+1];
                    geoArray[sp].psi  =angles[(i+dev)*3+2];
                    //precomute distances for faster execution
                    //Precompute per angle constant stuff for speed
                    computeDeltas(geoArray[sp],i+dev, &uvOrigin, &deltaU, &deltaV, &source);
                    
                    // compute upper bound limits of ray tracing
                    cropdist_init=maxdistanceCuboid(geo,i+dev); // TODO: this needs reworking for 3D

                    cudaSetDevice(dev);
                   
                   //TODO: we could do this around X and Y axis too, but we would need to compute the new axis of rotation (not possible to know from jsut the angles)
                    if (geo.theta==0.0f & geo.psi==0.0f){
                        kernelPixelDetector<false><<<grid,block,0,stream[dev*nStream_device]>>>(geoArray[sp],dProjection[(int)(i%(2*deviceCount)!=0)+dev*deviceCount], source, deltaU, deltaV, uvOrigin,geo.DSO[i+dev],floor(cropdist_init),texImg[dev]);
                    }
                    else{
                        kernelPixelDetector<true> <<<grid,block,0,stream[dev*nStream_device]>>>(geoArray[sp],dProjection[(int)(i%(2*deviceCount)!=0)+dev*deviceCount], source, deltaU, deltaV, uvOrigin,geo.DSO[i+dev],floor(cropdist_init),texImg[dev]);
                    }
                }
                
            }
            // Now that the computation is happening, we need to either prepare the memory for
            // combining of the projections (splits>1) or start removing previous results.
            
            
            // If our image does not fit in memory then we need to make sure we accumulate previous results too.
            // First, grab previous results and put them in the auxiliary variable
            if( !fits_in_memory && sp>0){
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    cudaMemcpyAsync(dProjection_accum[(int)(i%(2*deviceCount)!=0)+dev*deviceCount], result[i+dev], num_bytes_proj, cudaMemcpyHostToDevice,stream[dev*nStream_device+1]);                    
                }
            }
          
            // Second, take the results from previous compute call and add it to the code in execution.  
            if( !fits_in_memory && sp>0){
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    cudaStreamSynchronize(stream[dev*nStream_device+1]);
                    vecAddInPlaceInterp<<<(geo.nDetecU*geo.nDetecV+MAXTREADS-1)/MAXTREADS,MAXTREADS,0,stream[dev*nStream_device]>>>(dProjection[(int)(i%(2*deviceCount)!=0)+dev*deviceCount],dProjection_accum[(int)(i%(2*deviceCount)!=0)+dev*deviceCount],(unsigned long)geo.nDetecU*geo.nDetecV );
                    
                }
            }
            
 
            // Now, lets get out the projections from the previous execution of the kernels.
            if (i>0){
                for (dev = 0; dev < deviceCount; dev++){
                    
                    // copy result to host
                    cudaSetDevice(dev);
                    cudaMemcpyAsync(result[(i-deviceCount)+dev], dProjection[(int)(i%(2*deviceCount)==0)+dev*deviceCount], num_bytes_proj, cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
                }
            }
            
            // Make sure Computation on kernels has finished before we launch the next batch.
            for (dev = 0; dev < deviceCount; dev++){
                cudaSetDevice(dev);
                cudaStreamSynchronize(stream[dev*nStream_device]);
            }
            
        }
         
                // We still have the last one to get out, do that one
        
        int i = nangles-deviceCount+nangles%deviceCount;

        for (dev = 0; dev < deviceCount; dev++){
            if(i+dev<nangles){
                // copy result to host
                cudaSetDevice(dev);
                cudaDeviceSynchronize();
                cudaMemcpyAsync(result[i+dev], dProjection[(int)(i%(2*deviceCount)!=0)+dev*deviceCount], num_bytes_proj, cudaMemcpyDeviceToHost,stream[dev*nStream_device+1]);
            }
            
        }
        // Free memory for the next piece of image

        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(dev);
            cudaDestroyTextureObject(texImg[dev]);
            cudaFreeArray(d_cuArrTex[dev]);
            
        }
    }
    
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaFree(dProjection[dev*deviceCount]);
        cudaFree(dProjection[dev*deviceCount+1]);
        
    }
    free(dProjection);
    
    if(!fits_in_memory){
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(dev);
            cudaFree(dProjection_accum[dev*deviceCount]);
            cudaFree(dProjection_accum[dev*deviceCount+1]);
            
        }
        free(dProjection_accum);
    }
    freeGeoArray(splits,geoArray);

    cudaCheckErrors("cudaFree d_imagedata fail");
  
   for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]) ;

    cudaDeviceReset();
    return 0;
}
void CreateTextureInterp(int num_devices,const float* imagedata,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage)
{
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
     const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    for (unsigned int i = 0; i < num_devices; i++){
        cudaSetDevice(i);
        
        //cudaArray Descriptor
       
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        //cuda Array
        cudaMalloc3DArray(&d_cuArrTex[i], &channelDesc, extent);
        cudaCheckErrors("Texture memory allocation fail");
        
        
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
