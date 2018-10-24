/*-------------------------------------------------------------------------
 *
 * CUDA function for backrpojection using matched weigts for CBCT
 *
 *
 * CODE by  Ander Biguri & Sepideh Hatamikia
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
#include "voxel_backprojection2_spherical.hpp"
#include <stdio.h>
#include <math.h>

// https://stackoverflow.com/questions/16282136/is-there-a-cuda-equivalent-of-perror
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                printf("%s \n",msg);\
                printf("CBCT:CUDA:Atb",cudaGetErrorString(__err));\
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

__global__ void matrixConstantMultiply_sp(const Geometry geo,float* image,float constant){
    unsigned long long idx = threadIdx.x + blockIdx.x * blockDim.x;
     for(; idx<geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ; idx+=gridDim.x*blockDim.x) {
            image[idx]*=constant;
     }
    
}

// Using Matched weigths
__global__ void kernelPixelBackprojectionSpherical(const Geometry geo,
        float* image,
        int indAlpha,
        const float COR,
        const float DSD,
        const float DSO,
        const float cosa,
        const float sina,
        const Point3D deltaX ,
        const Point3D deltaY,
        const Point3D deltaZ,
        const Point3D xyzOrigin,           
        const Point3D xyzOffset,            // this is a direct copy, it has not been scaled
        const Point3D uv0Offset,
        Point3D source){           // This is a direct copy, it has not been scaled
    
    
    unsigned long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long indX = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long indZ = blockIdx.z * blockDim.z + threadIdx.z;
    //Make sure we dont go out of bounds
    unsigned long long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
    if (indX>=geo.nVoxelX | indY>=geo.nVoxelY |indZ>=geo.nVoxelZ)
        return;
    // Geometric trasnformations:
    
    //Source, scaled XYZ coordinates
   
    // "XYZ" in the scaled coordinate system of the current point. The iamge is rotated with the projection angles.
    Point3D P;
    P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
    P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y)-COR/geo.dDetecU;
    P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
        
    // This is the vector defining the line from the source to the Voxel
    float vectX,vectY,vectZ;
    vectX=(P.x -source.x);
    vectY=(P.y -source.y);
    vectZ=(P.z -source.z);
    
    
    // Get the coordinates in the detector UV where the mid point of the voxel is projected.
    float t=(DSO-DSD /*-DDO*/ - source.x)/vectX;
    float y,z;
    y=vectY*t+source.y;
    z=vectZ*t+source.z;
    float u,v;
    u=y+geo.nDetecU/2;
    v=z+geo.nDetecV/2;
    
    // TODO: put this in a separate kernel?
    // Compute the weigth of the matched backprojection , as in doi: 10.1088/0031-9155/56/13/004, eq (3)
    float weigth;
    //Real coordinates of Voxel. Instead of reverting the tranformation, its less math (faster) to compute it from the indexes.
    Point3D realvoxel; 
    realvoxel.x=-geo.sVoxelX/2+geo.dVoxelX/2    +indX*geo.dVoxelX   +xyzOffset.x;      
    realvoxel.y=-geo.sVoxelY/2+geo.dVoxelY/2    +indY*geo.dVoxelY   +xyzOffset.y; 
    realvoxel.z=-geo.sVoxelZ/2+geo.dVoxelZ/2    +indZ*geo.dVoxelZ   +xyzOffset.z;
    
    
    
    // This is not wrong. All numbers being computed here are for distances, (l,L), therefore
    // the entire rotation is not needed. I'd even say that the "alpha" is not needed either)
    
    
    //Real coords of Source
    // We already have S.x, and S.y and S.z are always zero. we just need to rotate
    Point3D S;
    S.x= DSO*cosa;
    S.y=-DSO*sina; 
   
    // Real XYZ coordinates of Detector.
    Point3D realD, realDaux; 
    // We know the index of the detector (u,v). Start from there.
    realDaux.x=-(DSD-DSO); 
    realDaux.y=-geo.sDetecU/2+geo.dDetecU/2 + u*geo.dDetecU +uv0Offset.x;
    realD.z   =-geo.sDetecV/2+geo.dDetecV/2 + v*geo.dDetecV +uv0Offset.y;
    //rotate the detector
    realD.x= realDaux.x*cosa  + realDaux.y*sina; //sin(-x)=-sin(x) , cos(-x)=cos(x)
    realD.y=-realDaux.x*sina  + realDaux.y*cosa; //sin(-x)=-sin(x) , cos(-x)=cos(x)
    float L,l;
    L = sqrt( (S.x-realD.x)*(S.x-realD.x)+ (S.y-realD.y)*(S.y-realD.y)+ (realD.z)*(realD.z)); // Sz=0 always.
    l = sqrt( (S.x-realvoxel.x)*(S.x-realvoxel.x)+ (S.y-realvoxel.y)*(S.y-realvoxel.y)+ (S.z-realvoxel.z)*(S.z-realvoxel.z));
    weigth=L*L*L/(DSD*l*l);
    
   // Get Value in the computed (U,V) and multiply by the corresponding weigth.
    image[idx]+=tex3D(tex, v  ,
                           u  ,
                           indAlpha+0.5)
                           *weigth;
    
    
}


int voxel_backprojection2_spherical(float const * const projections, Geometry geo, float* result,float const * const angles,int nalpha){
    
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
    
    
    Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, offDetec, source;
    
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
    for (unsigned int i=0;i<nalpha;i++){
        geo.alpha=-angles[i*3];
        geo.theta=-angles[i*3+1];
        geo.psi  =-angles[i*3+2];
        
        computeDeltasCubeSpherical(geo,i,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
        
        offOrig.x=geo.offOrigX[i];
        offOrig.y=geo.offOrigY[i];
        offOrig.z=geo.offOrigZ[i];
        offDetec.x=geo.offDetecU[i];
        offDetec.y=geo.offDetecV[i];
        
        kernelPixelBackprojectionSpherical<<<grid,block>>>
                (geo,dimage,i,geo.COR[i],geo.DSD[i],geo.DSO[i],cos(geo.alpha),sin(geo.alpha),deltaX,deltaY,deltaZ,xyzOrigin,offOrig,offDetec,source);
        cudaCheckErrors("Kernel fail");
    }
    // If we are timing this
    if (timekernel){
        cudaEventCreate(&stop);
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsedTime, start,stop);
        printf("%f\n" ,elapsedTime);
    }
     //in a Tesla, maximum blocks =15 SM * 4 blocks/SM=60
     matrixConstantMultiply_sp<<<60,MAXTREADS>>>( geo,dimage,geo.dVoxelX*geo.dVoxelY*geo.dVoxelZ/(geo.dDetecU*geo.dDetecV));
    
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
#ifndef BACKPROJECTION_SPHERICAL_HPP
void computeDeltasCubeSpherical(Geometry geo,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S){
    
    Point3D P, Px,Py,Pz;
    // Get coords of Img(0,0,0)
    P.x=-(geo.sVoxelX/2-geo.dVoxelX/2)+geo.offOrigX[i];
    P.y=-(geo.sVoxelY/2-geo.dVoxelY/2)+geo.offOrigY[i];
    P.z=-(geo.sVoxelZ/2-geo.dVoxelZ/2)+geo.offOrigZ[i];
    
    // Get coors from next voxel in each direction
    Px.x=P.x+geo.dVoxelX;      Py.x=P.x;                Pz.x=P.x;
    Px.y=P.y;                   Py.y=P.y+geo.dVoxelY;    Pz.y=P.y;
    Px.z=P.z;                   Py.z=P.z;                Pz.z=P.z+geo.dVoxelZ;
    
    
    
// Rotate image around X axis (this is equivalent of rotating the source and detector) RZ RY RZ
    
    eulerZYZ_back(geo,&P);
    eulerZYZ_back(geo,&Px);
    eulerZYZ_back(geo,&Py);
    eulerZYZ_back(geo,&Pz);

    
    
    //detector offset
    P.z =P.z-geo.offDetecV[i];            P.y =P.y-geo.offDetecU[i];
    Px.z =Px.z-geo.offDetecV[i];          Px.y =Px.y-geo.offDetecU[i];
    Py.z =Py.z-geo.offDetecV[i];          Py.y =Py.y-geo.offDetecU[i];
    Pz.z =Pz.z-geo.offDetecV[i];          Pz.y =Pz.y-geo.offDetecU[i];
    
    //Detector Roll pitch Yaw
    //
    //
    // first, we need to offset everything so (0,0,0) is the center of the detector
    // Only X is required for that
    P.x=P.x+(geo.DSD[i]-geo.DSO[i]);
    Px.x=Px.x+(geo.DSD[i]-geo.DSO[i]);
    Py.x=Py.x+(geo.DSD[i]-geo.DSO[i]);
    Pz.x=Pz.x+(geo.DSD[i]-geo.DSO[i]);
    rollPitchYawT(geo,i,&P);
    rollPitchYawT(geo,i,&Px);
    rollPitchYawT(geo,i,&Py);
    rollPitchYawT(geo,i,&Pz);
    
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);
    Px.x=Px.x-(geo.DSD[i]-geo.DSO[i]);
    Py.x=Py.x-(geo.DSD[i]-geo.DSO[i]);
    Pz.x=Pz.x-(geo.DSD[i]-geo.DSO[i]);
    //Done for P, now source
    
    source.x=geo.DSD[i]; //allready offseted for rotation
    source.y=-geo.offDetecU[i];
    source.z=-geo.offDetecV[i];
    rollPitchYawT(geo,i,&source);
    
    
    source.x=source.x-(geo.DSD[i]-geo.DSO[i]);//   source.y=source.y-auxOff.y;    source.z=source.z-auxOff.z;
    
//       mexPrintf("%f,%f,%f\n",source.x,source.y,source.z);
    // Scale coords so detector pixels are 1x1
    
    P.z =P.z /geo.dDetecV;                          P.y =P.y/geo.dDetecU;
    Px.z=Px.z/geo.dDetecV;                          Px.y=Px.y/geo.dDetecU;
    Py.z=Py.z/geo.dDetecV;                          Py.y=Py.y/geo.dDetecU;
    Pz.z=Pz.z/geo.dDetecV;                          Pz.y=Pz.y/geo.dDetecU;
    
    source.z=source.z/geo.dDetecV;                  source.y=source.y/geo.dDetecU;
    
    // get deltas of the changes in voxels
    deltaX->x=Px.x-P.x;   deltaX->y=Px.y-P.y;    deltaX->z=Px.z-P.z;
    deltaY->x=Py.x-P.x;   deltaY->y=Py.y-P.y;    deltaY->z=Py.z-P.z;
    deltaZ->x=Pz.x-P.x;   deltaZ->y=Pz.y-P.y;    deltaZ->z=Pz.z-P.z;
    
    
    *xyzorigin=P;
    *S=source;
}

#endif
