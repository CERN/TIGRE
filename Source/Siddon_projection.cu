/*
 * Code that uses texture memory to compute a 3D projection of CBCT
 *
 * IMPORTANT!!! CAUTION!! This code is designed for a Tesla 40k GPU.
 * It is a safe assumption to say that this code wont work in other GPUs as expected
 * or at all. Some of the involved reasons: float/double arithmetic.
 *
 * Ander Biguri
 */

#include <algorithm>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "Siddon_projection.hpp"
#include "mex.h"
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("CBCT:CUDA:Atb",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
// Declare the texture reference.
    texture<float, cudaTextureType3D , cudaReadModeElementType> tex;
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



__global__ void kernelPixelDetector( Geometry geo,
        float* detector,
        Point3D source ,
        Point3D deltaU,
        Point3D deltaV,
        Point3D uvOrigin){
    
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nDetecU* geo.nDetecV)
        return;
    
    
    
    
    /////// Get coordinates XYZ of pixel UV
    int pixelV = geo.nDetecV-idx % geo.nDetecV-1;
    int pixelU = idx / geo.nDetecV;
    Point3D pixel1D;
    pixel1D.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    pixel1D.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    pixel1D.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
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
    axm=min(-source.x/ray.x,(geo.nVoxelX-1-source.x)/ray.x);
    aym=min(-source.y/ray.y,(geo.nVoxelY-1-source.y)/ray.y);
    azm=min(-source.z/ray.z,(geo.nVoxelZ-1-source.z)/ray.z);
    axM=max(-source.x/ray.x,(geo.nVoxelX-1-source.x)/ray.x);
    ayM=max(-source.y/ray.y,(geo.nVoxelY-1-source.y)/ray.y);
    azM=max(-source.z/ray.z,(geo.nVoxelZ-1-source.z)/ray.z);
    float am=max(max(axm,aym),azm);
    float aM=min(min(axM,ayM),azM);
    
    // line intersects voxel space ->   am<aM
    if (am>=aM)
        detector[idx]=0;
    
    // Compute max/min image INDEX for intersection eq(11-19)
    // Discussion about ternary operator in CUDA: https://stackoverflow.com/questions/7104384/in-cuda-why-is-a-b010-more-efficient-than-an-if-else-version
    float imin,imax,jmin,jmax,kmin,kmax;
    // for X
    if( source.x<pixel1D.x){
        imin=(am==axm)? 1             : ceil (source.x+am*ray.x);
        imax=(aM==axM)? geo.nVoxelX-1 : floor(source.x+aM*ray.x);
    }else{
        imax=(am==axm)? geo.nVoxelX-2 : floor(source.x+am*ray.x);
        imin=(aM==axM)? 0             : ceil (source.x+aM*ray.x);
    }
    // for Y
    if( source.y<pixel1D.y){
        jmin=(am==aym)? 1             : ceil (source.y+am*ray.y);
        jmax=(aM==ayM)? geo.nVoxelY-1 : floor(source.y+aM*ray.y);
    }else{
        jmax=(am==aym)? geo.nVoxelY-2 : floor(source.y+am*ray.y);
        jmin=(aM==ayM)? 0             : ceil (source.y+aM*ray.y);
    }
    // for Z
    if( source.z<pixel1D.z){
        kmin=(am==azm)? 1             : ceil (source.z+am*ray.z);
        kmax=(aM==azM)? geo.nVoxelZ-1 : floor(source.z+aM*ray.z);
    }else{
        kmax=(am==azm)? geo.nVoxelZ-2 : floor(source.z+am*ray.z);
        kmin=(aM==azM)? 0             : ceil (source.z+aM*ray.z);
    }
    
    // get intersection point N1. eq(20-21) [(also eq 9-10)]
    float ax,ay,az;
    ax=(source.x<pixel1D.x)?  (imin-source.x)/ray.x  :  (imax-source.x)/ray.x;
    ay=(source.y<pixel1D.y)?  (jmin-source.y)/ray.y  :  (jmax-source.y)/ray.y;
    az=(source.z<pixel1D.z)?  (kmin-source.z)/ray.z  :  (kmax-source.z)/ray.z;
    
    
    
    // get index of first intersection. eq (26) and (19)
    int i,j,k;
    float aminc=min(min(ax,ay),az);
    i=(int)floor(source.x+ (aminc+am)/2*ray.x);
    j=(int)floor(source.y+ (aminc+am)/2*ray.y);
    k=(int)floor(source.z+ (aminc+am)/2*ray.z);
    // Initialize
    float ac=am;
    //eq (28), unit alphas
    float axu,ayu,azu;
    axu=1/abs(ray.x);
    ayu=1/abs(ray.y);
    azu=1/abs(ray.z);
    // eq(29), direction of update
    float iu,ju,ku;
    iu=(source.x< pixel1D.x)? 1 : -1;
    ju=(source.y< pixel1D.y)? 1 : -1;
    ku=(source.z< pixel1D.z)? 1 : -1;
    
    float maxlength=sqrt(ray.x*ray.x*geo.dVoxelX*geo.dVoxelX+ray.y*ray.y*geo.dVoxelY*geo.dVoxelY+ray.z*ray.z*geo.dVoxelZ*geo.dVoxelZ);
    float sum=0;
    int Np=(imax-imin+1)+(jmax-jmin+1)+(kmax-kmin+1); // Number of intersections
    // Go iterating over the line, intersection by intersection. If double point, no worries, 0 will be computed
//         while (ac<=aM){
    for (int ii=0;ii<Np;ii++){
        if (ax==aminc){
            sum+=(ax-ac)*tex3D(tex, i+0.5, j+0.5, k+0.5);
            i=i+iu;
            ac=ax;
            ax+=axu;
        }else if(ay==aminc){
            sum+=(ay-ac)*tex3D(tex, i+0.5, j+0.5, k+0.5);
//                 sum++;
            j=j+ju;
            ac=ay;
            ay+=ayu;
        }else if(az==aminc){
            sum+=(az-ac)*tex3D(tex, i+0.5, j+0.5, k+0.5);
//                 sum++;
            k=k+ku;
            ac=az;
            az+=azu;
        }
        aminc=min(min(ax,ay),az);
    }
    detector[idx]=sum*maxlength;
//         detector[idx]=source.y;
}


int siddon_ray_projection(float const * const img, Geometry geo, float** result,float const * const alphas,int nalpha){
    
    
    // BEFORE DOING ANYTHING: Use the proper CUDA enabled GPU: Tesla K40c
    
    // If you have another GPU and want to use this code, please change it, but make sure you know that is compatible.
    // also change MAXTREADS
    
//         int deviceCount = 0;
//         cudaGetDeviceCount(&deviceCount);
//         if (deviceCount == 0)
//         {
//             mexErrMsgIdAndTxt("CBCT:CUDA:Ax:cudaGetDeviceCount","No CUDA enabled NVIDIA GPUs found");
//         }
//         bool found=false;
//         for (int dev = 0; dev < deviceCount; ++dev)
//         {
//             cudaDeviceProp deviceProp;
//             cudaGetDeviceProperties(&deviceProp, dev);
//
//             if (strcmp(deviceProp.name, "Tesla K40c") == 0|| strcmp(deviceProp.name, "GeForce GT 740M") == 0){
//                 cudaSetDevice(dev);
//                 found=true;
//                 break;
//             }
//         }
//         if (!found)
//            mexErrMsgIdAndTxt("CBCT:CUDA:Ax:cudaDevice","No Supported GPU found");
    // DONE, Tesla found
    
    // copy data to CUDA memory
    cudaArray *d_imagedata = 0;
    
    const cudaExtent extent = make_cudaExtent(geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaMalloc3DArray(&d_imagedata, &channelDesc, extent);
    cudaCheckErrors("cudaMalloc3D error 3D tex");
    
    cudaMemcpy3DParms copyParams = { 0 };
    copyParams.srcPtr = make_cudaPitchedPtr((void*)img, extent.width*sizeof(float), extent.width, extent.height);
    copyParams.dstArray = d_imagedata;
    copyParams.extent = extent;
    copyParams.kind = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);
    
    cudaCheckErrors("cudaMemcpy3D fail");
    
    // Configure texture options
    tex.normalized = false;
    tex.filterMode = cudaFilterModePoint; //we dotn want itnerpolation
    tex.addressMode[0] = cudaAddressModeBorder;
    tex.addressMode[1] = cudaAddressModeBorder;
    tex.addressMode[2] = cudaAddressModeBorder;
    
    cudaBindTextureToArray(tex, d_imagedata, channelDesc);
    
    cudaCheckErrors("3D texture memory bind fail");
    
    
    
    
    //Done! Image put into texture memory.
    
    
    size_t num_bytes = geo.nDetecU*geo.nDetecV * sizeof(float);
    float* dProjection;
    cudaMalloc((void**)&dProjection, num_bytes);
    cudaCheckErrors("cudaMalloc fail");
    
    bool timekernel=false;
    cudaEvent_t start, stop;
    float elapsedTime;
    if (timekernel){
        cudaEventCreate(&start);
        cudaEventRecord(start,0);
    }
    Point3D source, deltaU, deltaV, uvOrigin;
    for (int i=0;i<nalpha;i++){
        
        geo.alpha=alphas[i];
        //precomute distances for faster execution
        //Precompute per angle constant stuff for speed
        computeDeltas_Siddon(geo,geo.alpha,i, &uvOrigin, &deltaU, &deltaV, &source);
        //Ray tracing!
        kernelPixelDetector<<<(geo.nDetecU*geo.nDetecV + MAXTREADS-1) / MAXTREADS,MAXTREADS>>>(geo,dProjection, source, deltaU, deltaV, uvOrigin);
        cudaCheckErrors("Kernel fail");
        // copy result to host
        cudaMemcpy(result[i], dProjection, num_bytes, cudaMemcpyDeviceToHost);
        cudaCheckErrors("cudaMemcpy fail");
        
        
    }
    if (timekernel){
        cudaEventCreate(&stop);
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsedTime, start,stop);
        mexPrintf("%f\n" ,elapsedTime);
    }
    
    cudaUnbindTexture(tex);
    cudaCheckErrors("Unbind  fail");
    cudaFree(dProjection);
    cudaFreeArray(d_imagedata);
    cudaCheckErrors("cudaFree d_imagedata fail");
    
    
    
    // tehre is no need to reset the device, but if one whants to use the NVIDIA Visual profiler, one should.
    cudaDeviceReset();
    return 0;
}



/* This code precomputes The location of the source and the Delta U and delta V (in the warped space)
 * to compute the locations of the x-rays. While it seems verbose and overly-optimized,
 * it does saves about 30% of each of the kernel calls. Thats something!
 **/
void computeDeltas_Siddon(Geometry geo, float alpha,int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){
    Point3D S;
    S.x=geo.DSO;
    S.y=0;
    S.z=0;
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD-geo.DSO);   P.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       P.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pu0.x=-(geo.DSD-geo.DSO);   Pu0.y= geo.dDetecU*(1-((float)geo.nDetecU/2)+0.5);       Pu0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pv0.x=-(geo.DSD-geo.DSO);   Pv0.y= geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       Pv0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-1);
    // Geomtric trasnformations:
    
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
    S2.x=S2.x-geo.offOrigX[i];       S2.y=S2.y-geo.offOrigY[i];       S2.z=S2.z-geo.offOrigZ[i];
    
    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2-geo.dVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2-geo.dVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2-geo.dVoxelZ/2;
    Pfinalu0.x=Pfinalu0.x+geo.sVoxelX/2-geo.dVoxelX/2;    Pfinalu0.y=Pfinalu0.y+geo.sVoxelY/2-geo.dVoxelY/2;        Pfinalu0.z=Pfinalu0.z+geo.sVoxelZ/2-geo.dVoxelZ/2;
    Pfinalv0.x=Pfinalv0.x+geo.sVoxelX/2-geo.dVoxelX/2;    Pfinalv0.y=Pfinalv0.y+geo.sVoxelY/2-geo.dVoxelY/2;        Pfinalv0.z=Pfinalv0.z+geo.sVoxelZ/2-geo.dVoxelZ/2;
    S2.x      =S2.x+geo.sVoxelX/2-geo.dVoxelX/2;          S2.y      =S2.y+geo.sVoxelY/2-geo.dVoxelY/2;              S2.z      =S2.z      +geo.sVoxelZ/2-geo.dVoxelZ/2;
    
    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S2.x      =S2.x/geo.dVoxelX;          S2.y      =S2.y/geo.dVoxelY;            S2.z      =S2.z/geo.dVoxelZ;
    
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
#ifndef PROJECTION_HPP

float maxDistanceCubeXY(Geometry geo, float alpha,int i){
    ///////////
    // Compute initial "t" so we access safely as less as out of bounds as possible.
    //////////
    
    
    float maxCubX,maxCubY;
    // Forgetting Z, compute max distance: diagonal+offset
    maxCubX=(geo.sVoxelX/2+ abs(geo.offOrigX[i]))/geo.dVoxelX;
    maxCubY=(geo.sVoxelY/2+ abs(geo.offOrigY[i]))/geo.dVoxelY;
    
    return geo.DSO/geo.dVoxelX-sqrt(maxCubX*maxCubX+maxCubY*maxCubY);
    
}
#endif
