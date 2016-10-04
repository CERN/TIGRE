/*-------------------------------------------------------------------------
 *
 * CUDA function for backrpojection using matched weigts for CBCT
 *
 *
 * CODE by  Ander Biguri
---------------------------------------------------------------------------
---------------------------------------------------------------------------
Copyright (c) 2015, University of Bath and CERN- European Organization for 
Nuclear Research
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
 ---------------------------------------------------------------------------

Contact: tigre.toolbox@gmail.com
Codes  : https://github.com/CERN/TIGRE
--------------------------------------------------------------------------- 
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
     *            |     /        /|             |
     *   A Z      |    /        / |*D           |
     *   |        |   +--------+  |             |
     *   |        |   |        |  |             |
     *   |        |   |     *O |  +             |
     *   *--->y   |   |        | /              |
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
    texture<float, cudaTextureType3D , cudaReadModeElementType> tex;

__global__ void matrixConstantMultiply(const Geometry geo,float* image,float constant){
    unsigned long long idx = threadIdx.x + blockIdx.x * blockDim.x;
     for(; idx<geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ; idx+=gridDim.x*blockDim.x) {
            image[idx]*=constant;
     }
    
}

// Using Matched weigths
__global__ void kernelPixelBackprojection(const Geometry geo,
        float* image,
        const int indAlpha,
        const Point3D deltaX ,
        const Point3D deltaY,
        const Point3D deltaZ,
        const Point3D xyzOrigin,           
        const Point3D xyzOffset,            // this is a direct copy, it has not been scaled
        const Point3D uv0Offset,           // This is a direct copy, it has not been scaled
        const float sinalpha,
        const float cosalpha){
    
    unsigned long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long indX = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long indZ = blockIdx.z * blockDim.z + threadIdx.z;
    //Make sure we dont go out of bounds
    unsigned long long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
    if (indX>=geo.nVoxelX | indY>=geo.nVoxelY |indZ>=geo.nVoxelZ)
        return;
    // Geometric trasnformations:
    
    //Source, scaled XYZ coordinates
    Point3D S;
    S.x=geo.DSO;                  // we dont scale the x direction, because the detecros is only in YZ (and the image is rotated)
    S.y=-uv0Offset.x/geo.dDetecU;            
    S.z=-uv0Offset.y/geo.dDetecV;
    // "XYZ" in the scaled coordinate system of the current point. The iamge is rotated with the projection angles.
    Point3D P;
    P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
    P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y)-geo.COR/geo.dDetecU;
    P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
        
    // This is the vector defining the line from the source to the Voxel
    float vectX,vectY,vectZ;
    vectX=(P.x -S.x);
    vectY=(P.y -S.y);
    vectZ=(P.z -S.z);
    
    
    // Get the coordinates in the detector UV where the mid point of the voxel is projected.
    float t=(geo.DSO-geo.DSD /*-DDO*/ - S.x)/vectX;
    float y,z;
    y=vectY*t+S.y;
    z=vectZ*t+S.z;
    float u,v;
    u=y+geo.nDetecU/2-0.5;
    v=z+geo.nDetecV/2-0.5;
    
    // TODO: put this in a separate kernel?
    // Compute the weigth of the matched backprojection , as in doi: 10.1088/0031-9155/56/13/004, eq (3)
    float weigth;
    //Real coordinates of Voxel. Instead of reverting the tranformation, its less math (faster) to compute it from the indexes.
    Point3D realvoxel; 
    realvoxel.x=-geo.sVoxelX/2+geo.dVoxelX/2    +indX*geo.dVoxelX   +xyzOffset.x;      
    realvoxel.y=-geo.sVoxelY/2+geo.dVoxelY/2    +indY*geo.dVoxelY   +xyzOffset.y; 
    realvoxel.z=-geo.sVoxelZ/2+geo.dVoxelZ/2    +indZ*geo.dVoxelZ   +xyzOffset.z;
    //Real coords of Source
    // We already have S.x, and S.y and S.z are always zero. we just need to rotate
    S.x= geo.DSO*cosalpha;
    S.y=-geo.DSO*sinalpha; 
   
    // Real XYZ coordinates of Detector.
    Point3D realD, realDaux; 
    // We know the index of the detector (u,v). Start from there.
    realDaux.x=-(geo.DSD-geo.DSO); 
    realDaux.y=-geo.sDetecU/2+geo.dDetecU/2 + u*geo.dDetecU +uv0Offset.x;
    realD.z   =-geo.sDetecV/2+geo.dDetecV/2 + v*geo.dDetecV +uv0Offset.y;
    //rotate the detector
    realD.x= realDaux.x*cosalpha  + realDaux.y*sinalpha; //sin(-x)=-sin(x) , cos(-x)=cos(x)
    realD.y=-realDaux.x*sinalpha  + realDaux.y*cosalpha; //sin(-x)=-sin(x) , cos(-x)=cos(x)
    float L,l;
    L = sqrt( (S.x-realD.x)*(S.x-realD.x)+ (S.y-realD.y)*(S.y-realD.y)+ (realD.z)*(realD.z)); // Sz=0 always.
    l = sqrt( (S.x-realvoxel.x)*(S.x-realvoxel.x)+ (S.y-realvoxel.y)*(S.y-realvoxel.y)+ (S.z-realvoxel.z)*(S.z-realvoxel.z));
    weigth=L*L*L/(geo.DSD*l*l);
    
   // Get Value in the computed (U,V) and multiply by the corresponding weigth.
    image[idx]+=tex3D(tex, u +0.5 ,
                           v +0.5 ,
                           indAlpha+0.5)
                           *weigth;
    
    
}


int voxel_backprojection2(float const * const projections, Geometry geo, float* result,float const * const alphas,int nalpha){
    
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
      int divx,divy,divz;
    
    divx=10;
    divy=10;
    divz=10;
    dim3 grid((geo.nVoxelX+divx-1)/divx,
              (geo.nVoxelY+divy-1)/divy,
              (geo.nVoxelZ+divz-1)/divz); 
    dim3 block(divx,divy,divz);
    // Main loop
    
    float sinalpha, cosalpha;
    for (unsigned int i=0;i<nalpha;i++){
        geo.alpha=-alphas[i];
        sinalpha=sin(geo.alpha);
        cosalpha=cos(geo.alpha);
        computeDeltasCube(geo,geo.alpha,i,&xyzOrigin,&deltaX,&deltaY,&deltaZ);
        
        offOrig.x=geo.offOrigX[i];
        offOrig.y=geo.offOrigY[i];
        offOrig.z=geo.offOrigZ[i];
        offDetec.x=geo.offDetecU[i];
        offDetec.y=geo.offDetecV[i];
        
        kernelPixelBackprojection<<<grid,block>>>
                (geo,dimage,i,deltaX,deltaY,deltaZ,xyzOrigin,offOrig,offDetec,sinalpha,cosalpha);
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
     //in a Tesla, maximum blocks =15 SM * 4 blocks/SM=60
     matrixConstantMultiply<<<60,MAXTREADS>>>( geo,dimage,geo.dVoxelX*geo.dVoxelY*geo.dVoxelZ/(geo.dDetecU*geo.dDetecV));
    
    cudaMemcpy(result, dimage, num_bytes, cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy fail");
    
    cudaUnbindTexture(tex);
    cudaCheckErrors("Unbind  fail");
    
    cudaFree(dimage);
    cudaFreeArray(d_projectiondata);
    cudaCheckErrors("cudaFree d_imagedata fail");
    //cudaDeviceReset();
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
    
    // Rotate image in the opposite direction of what the detector would rotate. We will keep detector still while the image
    // changes to accomodate any needed geometric transformation.
    
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
    
    // Compute unit vector of change between voxels
    deltaX->x=Px.x-P.x;   deltaX->y=Px.y-P.y;    deltaX->z=Px.z-P.z;
    deltaY->x=Py.x-P.x;   deltaY->y=Py.y-P.y;    deltaY->z=Py.z-P.z;
    deltaZ->x=Pz.x-P.x;   deltaZ->y=Pz.y-P.y;    deltaZ->z=Pz.z-P.z;
    
    // Detector offset is encoded in the image.
    P.z =P.z-geo.offDetecV[i]/geo.dDetecV;          
    P.y =P.y-geo.offDetecU[i]/geo.dDetecU;
    
    *xyzorigin=P;
    
}
#endif