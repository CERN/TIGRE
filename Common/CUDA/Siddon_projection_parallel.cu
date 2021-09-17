/*-------------------------------------------------------------------------
 *
 * CUDA functions for ray-voxel intersection based projection
 *
 * This file has the necesary fucntiosn to perform X-ray parallel projection
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
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("TIGRE:CUDA:Ax",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
// Declare the texture reference.
void CreateTextureParallel(float* image,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,cudaStream_t* stream);


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
    
    unsigned long u = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long v = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long projNumber=threadIdx.z;
            
    if (u>= geo.nDetecU || v>= geo.nDetecV || projNumber>=PROJ_PER_BLOCK)
        return;
    
    int indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    
#if IS_FOR_MATLAB_TIGRE
    size_t idx =  (size_t)(u  * geo.nDetecV + v)+ (size_t)projNumber*geo.nDetecV *geo.nDetecU ;
#else
    size_t idx =  (size_t)(v  * geo.nDetecU + u)+ (size_t)projNumber*geo.nDetecV *geo.nDetecU ;
#endif
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    

    /////// Get coordinates XYZ of pixel UV
    int pixelV = geo.nDetecV-v-1;
    int pixelU = u;
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
    // This variables are ommited because
    // bx,by,bz ={0,0,0}
    // dx,dy,dz ={1,1,1}
    // compute parameter values for x-ray parametric equation. eq(3-10)
    float axm,aym,azm;
    float axM,ayM,azM;
    
    /**************************************
     *
     *
     * Problem. In paralel beam, often ray.y or ray.x=0;
     * This leads to infinities progpagating and breaking everything.
     *
     * We need to fix it.
     *
     ***************************************/
    
    // In the paper Nx= number of X planes-> Nvoxel+1
    axm=fminf(-source.x/ray.x,(geo.nVoxelX-source.x)/ray.x);
    aym=fminf(-source.y/ray.y,(geo.nVoxelY-source.y)/ray.y);
//     azm=min(-source.z/ray.z,(geo.nVoxelZ-source.z)/ray.z);
    axM=fmaxf(-source.x/ray.x,(geo.nVoxelX-source.x)/ray.x);
    ayM=fmaxf(-source.y/ray.y,(geo.nVoxelY-source.y)/ray.y);
//     azM=max(-source.z/ray.z,(geo.nVoxelZ-source.z)/ray.z);
    float am=(fmaxf(axm,aym));
    float aM=(fminf(axM,ayM));
    
    // line intersects voxel space ->   am<aM
    if (am>=aM)
        detector[idx]=0.0f;
    
    // Compute max/min image INDEX for intersection eq(11-19)
    // Discussion about ternary operator in CUDA: https://stackoverflow.com/questions/7104384/in-cuda-why-is-a-b010-more-efficient-than-an-if-else-version
    float imin,imax,jmin,jmax;
    // for X
    if( source.x<pixel1D.x){
        imin=(am==axm)? 1.0f             : ceilf (source.x+am*ray.x);
        imax=(aM==axM)? geo.nVoxelX      : floorf(source.x+aM*ray.x);
    }else{
        imax=(am==axm)? geo.nVoxelX-1.0f : floorf(source.x+am*ray.x);
        imin=(aM==axM)? 0.0f             : ceilf (source.x+aM*ray.x);
    }
    // for Y
    if( source.y<pixel1D.y){
        jmin=(am==aym)? 1.0f             : ceilf (source.y+am*ray.y);
        jmax=(aM==ayM)? geo.nVoxelY      : floorf(source.y+aM*ray.y);
    }else{
        jmax=(am==aym)? geo.nVoxelY-1.0f : floorf(source.y+am*ray.y);
        jmin=(aM==ayM)? 0.0f             : ceilf (source.y+aM*ray.y);
    }
//     // for Z
//     if( source.z<pixel1D.z){
//         kmin=(am==azm)? 1             : ceilf (source.z+am*ray.z);
//         kmax=(aM==azM)? geo.nVoxelZ : floorf(source.z+aM*ray.z);
//     }else{
//         kmax=(am==azm)? geo.nVoxelZ-1 : floorf(source.z+am*ray.z);
//         kmin=(aM==azM)? 0             : ceilf (source.z+aM*ray.z);
//     }
    
    // get intersection point N1. eq(20-21) [(also eq 9-10)]
    float ax,ay;
    ax=(source.x<pixel1D.x)?  (imin-source.x)/(ray.x+0.000000000001f)  :  (imax-source.x)/(ray.x+0.000000000001f);
    ay=(source.y<pixel1D.y)?  (jmin-source.y)/(ray.y+0.000000000001f)  :  (jmax-source.y)/(ray.y+0.000000000001f);
//     az=(source.z<pixel1D.z)?  (kmin-source.z)/ray.z  :  (kmax-source.z)/ray.z;
    
    
    
    // get index of first intersection. eq (26) and (19)
    int i,j,k;
    float aminc=fminf(ax,ay);
    i=(int)floorf(source.x+ (aminc+am)/2*ray.x);
    j=(int)floorf(source.y+ (aminc+am)/2*ray.y);
    k=(int)floorf(source.z+ (aminc+am)/2*ray.z);
//     k=(int)source.z;
    // Initialize
    float ac=am;
    //eq (28), unit angles
    float axu,ayu;
    axu=1.0f/fabsf(ray.x);
    ayu=1.0f/fabsf(ray.y);
//     azu=1/abs(ray.z);
    // eq(29), direction of update
    float iu,ju;
    iu=(source.x< pixel1D.x)? 1.0f : -1.0f;
    ju=(source.y< pixel1D.y)? 1.0f : -1.0f;
//     ku=(source.z< pixel1D.z)? 1 : -1;
    
    float maxlength=sqrtf(ray.x*ray.x*geo.dVoxelX*geo.dVoxelX+ray.y*ray.y*geo.dVoxelY*geo.dVoxelY);//+ray.z*ray.z*geo.dVoxelZ*geo.dVoxelZ);
    float sum=0.0f;
    unsigned int Np=(imax-imin+1)+(jmax-jmin+1);//+(kmax-kmin+1); // Number of intersections
    // Go iterating over the line, intersection by intersection. If double point, no worries, 0 will be computed
    i+=0.5f;
    j+=0.5f;
    k+=0.5f;
    for (unsigned int ii=0;ii<Np;ii++){
        if (ax==aminc){
            sum+=(ax-ac)*tex3D<float>(tex, i, j, k);//(ax-ac)*
            i=i+iu;
            ac=ax;
            ax+=axu;
        }else if(ay==aminc){
            sum+=(ay-ac)*tex3D<float>(tex, i, j, k);//(ay-ac)*
            j=j+ju;
            ac=ay;
            ay+=ayu;
//         }else if(az==aminc){
//             sum+=(az-ac)*tex3D<float>(tex, i+0.5, j+0.5, k+0.5);
//             k=k+ku;
//             ac=az;
//             az+=azu;
        }
        aminc=fminf(ay,ax);
    }
    detector[idx]=maxlength*sum;
}


int siddon_ray_projection_parallel(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids){
    
    

    
    
    size_t num_bytes = (size_t)geo.nDetecU*(size_t)geo.nDetecV*(size_t)PROJ_PER_BLOCK* (size_t)sizeof(float);
    float** dProjection=(float **)malloc(2*sizeof(float *));
    for (int i = 0; i < 2; ++i){
        cudaMalloc((void**)&dProjection[i],   num_bytes);
        cudaCheckErrors("cudaMalloc projections fail");
    }
    int nStreams=2;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    
    for (int i = 0; i < 2; ++i){
        cudaStreamCreate(&stream[i]);
    }
    
    
        
    // Texture object variables
    cudaTextureObject_t *texImg = 0;
    cudaArray **d_cuArrTex = 0;
    texImg =(cudaTextureObject_t*)malloc(1*sizeof(cudaTextureObject_t));
    d_cuArrTex =(cudaArray**)malloc(1*sizeof(cudaArray*));
    
    CreateTextureParallel(img,geo,&d_cuArrTex[0], &texImg   [0],stream);
    cudaCheckErrors("Texture allocation fail");
    //Done! Image put into texture memory.

    
    
    Point3D source, deltaU, deltaV, uvOrigin;
    
    
    Point3D* projParamsArrayHost;
    cudaMallocHost((void**)&projParamsArrayHost,4*PROJ_PER_BLOCK*sizeof(Point3D));

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
            if(geo.alpha==0.0 || abs(geo.alpha-1.5707963267949)<0.0000001){
                geo.alpha=geo.alpha+1.1920929e-07;
            }
            
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
    
    
    //3: Rotate (around z)!
    Point3D Pfinal, Pfinalu0, Pfinalv0;
    
    Pfinal.x  =P.x*cos(geo.alpha)-P.y*sin(geo.alpha);       Pfinal.y  =P.y*cos(geo.alpha)+P.x*sin(geo.alpha);       Pfinal.z  =P.z;
    Pfinalu0.x=Pu0.x*cos(geo.alpha)-Pu0.y*sin(geo.alpha);   Pfinalu0.y=Pu0.y*cos(geo.alpha)+Pu0.x*sin(geo.alpha);   Pfinalu0.z=Pu0.z;
    Pfinalv0.x=Pv0.x*cos(geo.alpha)-Pv0.y*sin(geo.alpha);   Pfinalv0.y=Pv0.y*cos(geo.alpha)+Pv0.x*sin(geo.alpha);   Pfinalv0.z=Pv0.z;
    
    Point3D S2;
    S2.x=S.x*cos(geo.alpha)-S.y*sin(geo.alpha);
    S2.y=S.y*cos(geo.alpha)+S.x*sin(geo.alpha);
    S2.z=S.z;
    
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
void CreateTextureParallel(float* image,Geometry geo,cudaArray** d_cuArrTex, cudaTextureObject_t *texImage,cudaStream_t* stream){    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    
    
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
