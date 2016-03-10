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
#include "voxel_backprojection.hpp"
#include "voxel_backprojection2.hpp"
#include "mex.h"
#include <math.h>

// https://stackoverflow.com/questions/16282136/is-there-a-cuda-equivalent-of-perror
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("CBCT:CUDA:Atb",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
    
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
              |     /        /|             |
     A Z      |    /        / |*D           |
     |        |   +--------+  |             |
     |        |   |        |  |             |
     |        |   |     *O |  +             |
     *--->y   |   |        | /              |
    /         |   |        |/               |
   V X        |   +--------+                |
     *            |-----------------------------|
     *
     *           *S
     *
     *
     *
     *
     *
     **/
    texture<float, cudaTextureType3D , cudaReadModeElementType> tex;

__global__ void matrixConstantMultiply(const Geometry geo,float* image,float constant){
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ)
        return;
    image[idx]*=constant;
    
}
__global__ void kernelPixelBackprojection(const Geometry geo,
        float* image,
        int indAlpha,
        Point3D deltaX ,
        Point3D deltaY,
        Point3D deltaZ,
        Point3D xyzOrigin,
        Point3D xyzOffset,
        Point3D uv0Offset){
    //Make sure we dont go out of bounds
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ)
        return;
    
    
    int indZ = idx / (geo.nVoxelY * geo.nVoxelX);
    int resZ = idx % (geo.nVoxelY * geo.nVoxelX);
    int indY = resZ / geo.nVoxelX;
    int indX = resZ % geo.nVoxelX;
    // Geometric trasnformations:
    
    //Source
    Point3D S;
    S.x=geo.DSO;
    S.y=-uv0Offset.x;
    S.z=-uv0Offset.y/geo.dDetecV;
    // "XYZ" in the warped coordinate system of the current point
    Point3D P;
    P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
    P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y);
    P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
    
    
    // compute the weigth for the backprojection. This needs the X and Y coords on the real workd of the image
    
    
    float vectX,vectY,vectZ;
    vectX=(P.x -S.x);
    vectY=(P.y -S.y);
    vectZ=(P.z -S.z);
    
    
    // Get the coordinates in the projection where the mid point of the voxel is projected.
    float t=(geo.DSO-geo.DSD /*-DDO*/ - S.x)/vectX;
    float y,z;
    y=vectY*t+S.y;
    z=vectZ*t+S.z;
    float u,v;
    u=y+geo.nDetecU/2-0.5;
    v=z+geo.nDetecV/2-0.5;
    
    
    // Get interpolated value in the current projection
    float weigth;
    float realx,realy,realz; //of voxel
    realx=-geo.sVoxelX/2+geo.dVoxelX/2    +indX*geo.dVoxelX   -xyzOffset.x; // /geo.dDetecU;  X never gets scaled.
    realy=-geo.sVoxelY/2+geo.dVoxelY/2    +indY*geo.dVoxelY   -xyzOffset.y/geo.dDetecU; // and Y gets scalled by U
    realz=-geo.sVoxelZ/2+geo.dVoxelZ/2    +indZ*geo.dVoxelZ   -xyzOffset.z/geo.dDetecV;
    float realDx,realDy,realDz; //of detector
    realDx=-(geo.DSD-geo.DSO);
    realDy=y*geo.dDetecU;
    realDz=z*geo.dDetecV;
    float L,l;
    L =sqrt( (S.x-realDx)*(S.x-realDx)+ (S.y-realDy)*(S.y-realDy)+ (S.z-realDz)*(S.z-realDz));
    l=sqrt( (S.x-realx)*(S.x-realx)+ (S.y-realy)*(S.y-realy)+ (S.z-realz)*(S.z-realz));
    weigth=L*L*L/(geo.DSD*l*l);
    image[idx]+=tex3D(tex, u +0.5 ,
            v +0.5 ,
            indAlpha                                           +0.5)*weigth;
    //test this and implement the deltas
    
    
}


int voxel_backprojection2(float const * const projections, Geometry geo, float* result,float const * const alphas,int nalpha){
    
//     // BEFORE DOING ANYTHING: Use the proper CUDA enabled GPU: Tesla K40c
//     int deviceCount = 0;
//     cudaGetDeviceCount(&deviceCount);
//      if (deviceCount == 0)
//     {
//         mexErrMsgIdAndTxt("CBCT:CUDA:Atb:cudaGetDeviceCount","No CUDA enabled NVIDIA GPUs found");
//     }
//     bool found=false;
//     for (int dev = 0; dev < deviceCount; ++dev)
//     {
//         cudaDeviceProp deviceProp;
//         cudaGetDeviceProperties(&deviceProp, dev);
//
//         if (strcmp(deviceProp.name, "Tesla K40c") == 0 || strcmp(deviceProp.name, "GeForce GT 740M") == 0){
//             cudaSetDevice(dev);
//             found=true;
//             break;
//         }
//     }
//     if (!found)
//        mexErrMsgIdAndTxt("CBCT:CUDA:Ax:cudaDevice","No Supported GPU found");
    
    
    // Done, Tesla found.
    
    /*
     * Allocate texture memory on the device
     */
    
    // copy data to CUDA memory
    cudaArray *d_projectiondata = 0;
    const cudaExtent extent = make_cudaExtent(geo.nDetecU,geo.nDetecV,nalpha);
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaMalloc3DArray(&d_projectiondata, &channelDesc, extent);
    cudaCheckErrors("cudaMalloc3D error 3D tex");
    
    cudaMemcpy3DParms copyParams = { 0 };
    copyParams.srcPtr = make_cudaPitchedPtr((void*)projections, extent.width*sizeof(float), extent.width, extent.height);
    copyParams.dstArray = d_projectiondata;
    copyParams.extent = extent;
    copyParams.kind = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);
    
    cudaCheckErrors("cudaMemcpy3D fail");
    
    // Configure texture options
    tex.normalized = false;
    tex.filterMode = cudaFilterModeLinear;
    tex.addressMode[0] = cudaAddressModeBorder;
    tex.addressMode[1] = cudaAddressModeBorder;
    tex.addressMode[2] = cudaAddressModeBorder;
    
    cudaBindTextureToArray(tex, d_projectiondata, channelDesc);
    
    cudaCheckErrors("3D texture memory bind fail");
    
    
    // Allocate result image memory
    size_t num_bytes = geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ * sizeof(float);
    float* dimage;
    cudaMalloc((void**)&dimage, num_bytes);
    cudaMemset(dimage,0,num_bytes);
    cudaCheckErrors("cudaMalloc fail");
    
    
    Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, offDetec;
    
    // time the kernel?
    bool timekernel=false;
    cudaEvent_t start, stop;
    float elapsedTime;
    if (timekernel){
        
        cudaEventCreate(&start);
        cudaEventRecord(start,0);
    }
    
    for (int i=0;i<nalpha;i++){
        geo.alpha=-alphas[i];
        computeDeltasCube(geo,geo.alpha,i,&xyzOrigin,&deltaX,&deltaY,&deltaZ);
        
        offOrig.x=geo.offOrigX[i];
        offOrig.y=geo.offOrigY[i];
        offDetec.x=geo.offDetecU[i];
        offDetec.y=geo.offDetecV[i];
        
        kernelPixelBackprojection<<<(geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ + MAXTREADS-1) / MAXTREADS,MAXTREADS>>>
                (geo,dimage,i,deltaX,deltaY,deltaZ,xyzOrigin,offOrig,offDetec);
        cudaCheckErrors("Kernel fail");
    }
    // If we are timing this
    if (timekernel){
        cudaEventCreate(&stop);
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsedTime, start,stop);
        mexPrintf("%f\n" ,elapsedTime);
    }
    matrixConstantMultiply<<<(geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ + MAXTREADS-1) / MAXTREADS,MAXTREADS>>>( geo,dimage,geo.dVoxelX*geo.dVoxelY*geo.dVoxelZ/(geo.dDetecU*geo.dDetecV));
    
    cudaMemcpy(result, dimage, num_bytes, cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy fail");
    
    cudaUnbindTexture(tex);
    cudaCheckErrors("Unbind  fail");
    
    cudaFree(dimage);
    cudaFreeArray(d_projectiondata);
    cudaCheckErrors("cudaFree d_imagedata fail");
    cudaDeviceReset();
    return 0;
    
}
#ifndef BACKPROJECTION_HPP
void computeDeltasCube(Geometry geo, float alpha,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ){
    
    Point3D P0, Px0,Py0,Pz0;
    // Get coords of Img(0,0,0)
    P0.x=-(geo.sVoxelX/2-geo.dVoxelX/2)+geo.offOrigX[i];
    P0.y=-(geo.sVoxelY/2-geo.dVoxelY/2)+geo.offOrigY[i];
    P0.z=-(geo.sVoxelZ/2-geo.dVoxelZ/2)+geo.offOrigZ[i];
    
    // Get coors from next voxel in each direction
    Px0.x=P0.x+geo.dVoxelX;       Py0.x=P0.x;                Pz0.x=P0.x;
    Px0.y=P0.y;                   Py0.y=P0.y+geo.dVoxelY;    Pz0.y=P0.y;
    Px0.z=P0.z;                   Py0.z=P0.z;                Pz0.z=P0.z+geo.dVoxelZ;
    
    // Rotate image (this is equivalent of rotating the source and detector)
    
    Point3D P, Px,Py,Pz; // We need other auxiliar variables to be able to perform the rotation, or we would overwrite values!
    P.x =P0.x *cos(alpha)-P0.y *sin(alpha);       P.y =P0.x *sin(alpha)+P0.y *cos(alpha);      P.z =P0.z;
    Px.x=Px0.x*cos(alpha)-Px0.y*sin(alpha);       Px.y=Px0.x*sin(alpha)+Px0.y*cos(alpha);      Px.z=Px0.z;
    Py.x=Py0.x*cos(alpha)-Py0.y*sin(alpha);       Py.y=Py0.x*sin(alpha)+Py0.y*cos(alpha);      Py.z=Py0.z;
    Pz.x=Pz0.x*cos(alpha)-Pz0.y*sin(alpha);       Pz.y=Pz0.x*sin(alpha)+Pz0.y*cos(alpha);      Pz.z=Pz0.z;
    
    // Scale coords so detector pixels are 1x1
    
    P.z =P.z /geo.dDetecV;                          P.y =P.y/geo.dDetecU;
    Px.z=Px.z/geo.dDetecV;                          Px.y=Px.y/geo.dDetecU;
    Py.z=Py.z/geo.dDetecV;                          Py.y=Py.y/geo.dDetecU;
    Pz.z=Pz.z/geo.dDetecV;                          Pz.y=Pz.y/geo.dDetecU;
    
    
    deltaX->x=Px.x-P.x;   deltaX->y=Px.y-P.y;    deltaX->z=Px.z-P.z;
    deltaY->x=Py.x-P.x;   deltaY->y=Py.y-P.y;    deltaY->z=Py.z-P.z;
    deltaZ->x=Pz.x-P.x;   deltaZ->y=Pz.y-P.y;    deltaZ->z=Pz.z-P.z;
    
    
    P.z =P.z-geo.offDetecV[i]/geo.dDetecV;          P.y =P.y-geo.offDetecU[i]/geo.dDetecU;
    *xyzorigin=P;
    
}
#endif