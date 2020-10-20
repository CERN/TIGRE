/*-------------------------------------------------------------------------
 *
 * CUDA functions for texture-memory interpolation based projection
 *
 * This file has the necesary functions to perform X-ray parallel projection
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
#include "ray_interpolated_projection_parallel.hpp"
#include "mex.h"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("TIGRE:Ax:interpolated_parallel",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
// Declare the texture reference.
    texture<float, cudaTextureType3D , cudaReadModeElementType> tex;

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
void CreateTextureParallelInterp(float* image,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,cudaStream_t* stream);
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
    
    unsigned long y = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long x = blockIdx.x * blockDim.x + threadIdx.x;
//     unsigned long idx =  x  * geo.nDetecV + y;
    unsigned long projNumber=threadIdx.z;
    
    if ((x>= geo.nDetecU) | (y>= geo.nDetecV)|  (projNumber>=PROJ_PER_BLOCK))
        return;
    
    int indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    
    size_t idx =  (size_t)(x  * geo.nDetecV + y)+ (size_t)projNumber*geo.nDetecV *geo.nDetecU ;
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    
    float DSO = projFloatsArrayDev[2*projNumber+0];
    float maxdist = projFloatsArrayDev[2*projNumber+1];
    
    
    /////// Get coordinates XYZ of pixel UV
    int pixelV = geo.nDetecV-y-1;
    int pixelU = x;
    
    
    
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



int interpolation_projection_parallel(float  *  img, Geometry geo, float** result,float const * const angles,int nangles){
    
    
    
    size_t num_bytes = geo.nDetecU*geo.nDetecV *PROJ_PER_BLOCK* sizeof(float);
    float** dProjection=(float **)malloc(2*sizeof(float *));
    for (int i = 0; i < 2; ++i){
        cudaMalloc((void**)&dProjection[i],   num_bytes);
        cudaCheckErrors("cudaMalloc projections fail");
    }
    // allocate streams for memory and compute
    int nStreams=2;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));;
    
    for (int i = 0; i < 2; ++i){
        cudaStreamCreate(&stream[i]);
    }
    
    
    // Texture object variables
    cudaTextureObject_t *texImg = 0;
    cudaArray **d_cuArrTex = 0;
    texImg =(cudaTextureObject_t*)malloc(1*sizeof(cudaTextureObject_t));
    d_cuArrTex =(cudaArray**)malloc(1*sizeof(cudaArray*));
    
    CreateTextureParallelInterp(img,geo,&d_cuArrTex[0], &texImg[0],stream);
    cudaCheckErrors("Texture allocation fail");
    //Done! Image put into texture memory.
    
    
    
    Point3D source, deltaU, deltaV, uvOrigin;
    
    Point3D* projParamsArrayHost;
    cudaMallocHost((void**)&projParamsArrayHost,4*PROJ_PER_BLOCK*sizeof(Point3D));
    float* projFloatsArrayHost;
    cudaMallocHost((void**)&projFloatsArrayHost,2*PROJ_PER_BLOCK*sizeof(float));
    
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
    
    float maxdist;
    for ( i=0; i<noOfKernelCalls; i++){
        for(unsigned int j=0; j<PROJ_PER_BLOCK; j++){
            proj_global=i*PROJ_PER_BLOCK+j;
            if (proj_global>=nangles)
                break;
            
            geo.alpha=angles[proj_global*3];
            geo.theta=angles[proj_global*3+1];
            geo.psi  =angles[proj_global*3+2];
            //precomute distances for faster execution
            maxdist=maxdistanceCuboid(geo,proj_global);
            //Precompute per angle constant stuff for speed
            computeDeltas_parallel(geo,geo.alpha,proj_global, &uvOrigin, &deltaU, &deltaV, &source);
            //Ray tracing!
            projParamsArrayHost[4*j]=uvOrigin;		// 6*j because we have 6 Point3D values per projection
            projParamsArrayHost[4*j+1]=deltaU;
            projParamsArrayHost[4*j+2]=deltaV;
            projParamsArrayHost[4*j+3]=source;
            
            projFloatsArrayHost[2*j]=geo.DSO[proj_global];
            projFloatsArrayHost[2*j+1]=floor(maxdist);
            
        }
        cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*4*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[0]);
        cudaMemcpyToSymbolAsync(projFloatsArrayDev, projFloatsArrayHost, sizeof(float)*2*PROJ_PER_BLOCK,0,cudaMemcpyHostToDevice,stream[0]);
        cudaStreamSynchronize(stream[0]);
        
        kernelPixelDetector_parallel_interpolated<<<numBlocks,threadsPerBlock,0,stream[0]>>>(geo,dProjection[(int)i%2==0],i,nangles,texImg[0]);
        // copy result to host
        if (i>0)
             cudaMemcpyAsync(result[i*PROJ_PER_BLOCK-PROJ_PER_BLOCK],dProjection[(int)i%2!=0], num_bytes, cudaMemcpyDeviceToHost,stream[1]);    
    }
    cudaDeviceSynchronize();
    
    int lastangles=nangles-(i-1)*PROJ_PER_BLOCK;
    cudaMemcpyAsync(result[(i-1)*PROJ_PER_BLOCK],dProjection[(int)(i-1)%2==0], lastangles*geo.nDetecV*geo.nDetecU*sizeof(float), cudaMemcpyDeviceToHost,stream[1]);

    
    cudaDestroyTextureObject(texImg[0]);
    cudaFreeArray(d_cuArrTex[0]);
    free(texImg); texImg = 0;
    free(d_cuArrTex); d_cuArrTex = 0;
    cudaCheckErrors("Unbind  fail");
    cudaFree(dProjection[0]);
    cudaFree(dProjection[1]);
    free(dProjection);
    cudaFreeHost(projParamsArrayHost);
    cudaFreeHost(projFloatsArrayHost);

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
    // Geomtric trasnformations:
    
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
void CreateTextureParallelInterp(float* image,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,cudaStream_t* stream){    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    
    
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
    copyParams.kind     = cudaMemcpyHostToDevice;
    cudaMemcpy3DAsync(&copyParams,stream[1]);
    
    
    //Array creation End
    
    cudaResourceDesc    texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array  = d_cuArrTex[0];
    cudaTextureDesc     texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));
    texDescr.normalizedCoords = false;
    texDescr.filterMode = cudaFilterModeLinear;
    texDescr.addressMode[0] = cudaAddressModeBorder;
    texDescr.addressMode[1] = cudaAddressModeBorder;
    texDescr.addressMode[2] = cudaAddressModeBorder;
    texDescr.readMode = cudaReadModeElementType;
    cudaCreateTextureObject(&texImage[0], &texRes, &texDescr, NULL);
    
}