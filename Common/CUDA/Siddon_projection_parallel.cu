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
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "Siddon_projection_parallel.hpp"
#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                tigreCudaError("TIGRE:CUDA:Ax", msg, cudaGetErrorString(__err));\
        } \
} while (0)
    
    
// Declare the texture reference.
void CreateTextureParallel(float* image,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,cudaStream_t* stream, cudaMemcpyKind kind);


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
    detector[idx]=maxlength*sum;
}


int siddon_ray_projection_parallel(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // With GPU data the copies to the caller's array are device-to-device and
    // asynchronous, so they go on the compute stream: that orders them with the
    // kernels that reuse the double-buffered projections. (Copies to pageable
    // host memory are synchronous, so the copy stream is safe for them.)
    const int copyStream = gpuData ? 0 : 1;
    // This projector runs on the current device; GPU data must be used on its own.
    if (gpuData) cudaSetDevice(gpuids[0]);
    const cudaMemcpyKind kindOut = gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost;
    
    

    
    
    size_t num_bytes = (size_t)geo.nDetecU*(size_t)geo.nDetecV*(size_t)PROJ_PER_BLOCK* (size_t)sizeof(float);
    float** dProjection=(float **)malloc(2*sizeof(float *));
    dProjection[0] = 0; dProjection[1] = 0;
    // If a CUDA error raises (Python build) before the end of this function,
    // onError releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    onError.add([dProjection]{ cudaFree(dProjection[0]); cudaFree(dProjection[1]); free(dProjection); });
    for (int i = 0; i < 2; ++i){
        cudaMalloc((void**)&dProjection[i],   num_bytes);
        cudaCheckErrors("cudaMalloc projections fail");
    }
    int nStreams=2;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    memset(stream, 0, nStreams*sizeof(cudaStream_t));
    onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
    
    for (int i = 0; i < 2; ++i){
        cudaStreamCreate(&stream[i]);
    }
    
    
        
    // Texture object variables
    cudaTextureObject_t *texImg = 0;
    cudaArray **d_cuArrTex = 0;
    texImg =(cudaTextureObject_t*)malloc(1*sizeof(cudaTextureObject_t));
    d_cuArrTex =(cudaArray**)malloc(1*sizeof(cudaArray*));
    texImg[0] = 0; d_cuArrTex[0] = 0;
    onError.add([texImg, d_cuArrTex]{ if (texImg[0]) cudaDestroyTextureObject(texImg[0]); if (d_cuArrTex[0]) cudaFreeArray(d_cuArrTex[0]); free(texImg); free(d_cuArrTex); });
    
    CreateTextureParallel(img,geo,&d_cuArrTex[0], &texImg   [0],stream,(gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice));
    cudaCheckErrors("Texture allocation fail");
    //Done! Image put into texture memory.

    
    
    Point3D source, deltaU, deltaV, uvOrigin;
    
    
    Point3D* projParamsArrayHost = 0;
    cudaMallocHost((void**)&projParamsArrayHost,4*PROJ_PER_BLOCK*sizeof(Point3D));
    onError.add([projParamsArrayHost]{ cudaFreeHost(projParamsArrayHost); });

    // 16x16 gave the best performance empirically
    // Funnily that makes it compatible with most GPUs.....
    int divU,divV,divangle;
    divU=PIXEL_SIZE_BLOCK;
    divV=PIXEL_SIZE_BLOCK;
    
    dim3 numBlocks((geo.nDetecU+divU-1)/divU,(geo.nDetecV+divV-1)/divV,1);
    
    dim3 threadsPerBlock(divU,divV,PROJ_PER_BLOCK);
    
    unsigned int proj_global;
    unsigned int noOfKernelCalls = (nangles+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_BLOCK
    unsigned int i;
    for ( i=0; i<noOfKernelCalls; i++){
        
         for(unsigned int j=0; j<PROJ_PER_BLOCK; j++){
            proj_global=i*PROJ_PER_BLOCK+j;
            if (proj_global>=nangles)
               break;
            geo.alpha=angles[proj_global*3];
            geo.theta=angles[proj_global*3+1];
            geo.psi  =angles[proj_global*3+2];
            
            //precomute distances for faster execution
            //Precompute per angle constant stuff for speed
            computeDeltas_Siddon_parallel(geo,geo.alpha,proj_global, &uvOrigin, &deltaU, &deltaV, &source);
            //Ray tracing!
            projParamsArrayHost[4*j]=uvOrigin;		// 6*j because we have 6 Point3D values per projection
            projParamsArrayHost[4*j+1]=deltaU;
            projParamsArrayHost[4*j+2]=deltaV;
            projParamsArrayHost[4*j+3]=source;

         }
         
         cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[0]);
         cudaStreamSynchronize(stream[0]);
         kernelPixelDetector_parallel<<<numBlocks,threadsPerBlock,0,stream[0]>>>(geo,dProjection[(int)i%2==0],i,nangles,texImg[0]);
         // copy result to host
         if (i>0)
             cudaMemcpyAsync(result[i*PROJ_PER_BLOCK-PROJ_PER_BLOCK],dProjection[(int)i%2!=0], num_bytes, kindOut,stream[copyStream]);
    }
    cudaDeviceSynchronize();
    
    int lastangles=nangles-(i-1)*PROJ_PER_BLOCK;
    cudaMemcpyAsync(result[(i-1)*PROJ_PER_BLOCK],dProjection[(int)(i-1)%2==0], lastangles*geo.nDetecV*geo.nDetecU*sizeof(float), kindOut,stream[copyStream]);

    

    if (gpuData) cudaDeviceSynchronize();  // the caller reads the result as soon as we return
    cudaCheckErrors("Projection fail");
    onError.dismiss();  // from here on the function releases everything itself
    cudaDestroyTextureObject(texImg[0]);
    cudaFreeArray(d_cuArrTex[0]);
    free(texImg); texImg = 0;
    free(d_cuArrTex); d_cuArrTex = 0;
    cudaCheckErrors("Unbind  fail");
    cudaFree(dProjection[0]);
    cudaFree(dProjection[1]);
    free(dProjection);
    cudaFreeHost(projParamsArrayHost);
    cudaCheckErrors("cudaFree d_imagedata fail");
    
    
    for (int i = 0; i < 2; ++i){
      cudaStreamDestroy(stream[i]);
    }
//     cudaDeviceReset();
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
void CreateTextureParallel(float* image,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,cudaStream_t* stream, cudaMemcpyKind kind){    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    
    
    const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
  
    //cudaArray Descriptor
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    //cuda Array
    cudaMalloc3DArray(&d_cuArrTex[0], &channelDesc, extent);


        cudaMemcpy3DParms copyParams = {0};
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)image, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[0];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
        cudaMemcpy3DAsync(&copyParams,stream[1]);
        // Wait for the copy: the kernels run on stream[0] and would otherwise read
        // the texture while it is still being filled. This raced always for GPU
        // arrays (device-to-device copies are fully asynchronous), and for host
        // data whenever the copy outlived the call (pinned memory, or the last
        // DMA chunk of a pageable copy).
        cudaStreamSynchronize(stream[1]);
    

    //Array creation End

        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArrTex[0];
        cudaTextureDesc     texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = false;
        texDescr.filterMode = cudaFilterModePoint;
        texDescr.addressMode[0] = cudaAddressModeBorder;
        texDescr.addressMode[1] = cudaAddressModeBorder;
        texDescr.addressMode[2] = cudaAddressModeBorder;
        texDescr.readMode = cudaReadModeElementType;
        cudaCreateTextureObject(&texImage[0], &texRes, &texDescr, NULL);
    
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
