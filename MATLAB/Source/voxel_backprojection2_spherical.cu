/*-------------------------------------------------------------------------
 *
 * CUDA function for backrpojection using matched weigts for CBCT
 *
 *
 * CODE by  Ander Biguri & Sepideh Hatamikia
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
#include "voxel_backprojection.hpp"
#include "voxel_backprojection2.hpp"
#include "voxel_backprojection2_spherical.hpp"
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
    const int VOXELS_PER_THREAD = 1;
void CreateTexture2Spherical(int num_devices,const float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage);

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
        Point3D source,
        cudaTextureObject_t tex){           // This is a direct copy, it has not been scaled
    
    
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
    
    float sample=tex3D<float>(tex, v ,u , indAlpha+0.5f);
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
    image[idx]+=sample*weigth;
    
    
}


int voxel_backprojection2_spherical(float const * const projections, Geometry geo, float* result,float const * const angles,int nalpha){
    
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    cudaCheckErrors("Device query fail");
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Atb:Voxel_backprojection_spherical:GPUselect","There are no available device(s) that support CUDA\n");
    }
    
    // Check the available devices, and if they are the same
    int dev;
    checkDevices();
    
    // Split the CT problem
    unsigned int split_image;
    unsigned int split_projections;
    splitCTbackprojection(deviceCount,geo,nalpha,&split_image,&split_projections);
    
    // Create the arrays for the geometry. The main difference is that geo.offZ has been tuned for the
    // image slices. The rest of the Geometry is the same
    Geometry* geoArray=(Geometry*)malloc(split_image*deviceCount*sizeof(Geometry));
    createGeoArray(split_image*deviceCount,geo,geoArray,nalpha);
    
    // Now lest allocate all the image memory on the GPU, so we can use it later. If we have made our numbers correctly
    // in the previous section this should leave enough space for the textures.
    size_t num_bytes_img = geo.nVoxelX*geo.nVoxelY*geoArray[0].nVoxelZ* sizeof(float);
    
    
    float** dimage=(float**)malloc(deviceCount*sizeof(float*));
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaMalloc((void**)&dimage[dev], num_bytes_img);
        cudaCheckErrors("cudaMalloc fail");
    }
    
    // Start with the main loop. The Projection data needs to be allocated and dealocated in the main loop
    // as due to the nature of cudaArrays, we can not reuse them. This should not be a problem for the fast execution
    // of the code, as repeated allocation and deallocation only happens when the projection data is very very big,
    // and therefore allcoation time should be negligible, fluctuation of other computations should mask the time.
    unsigned long long proj_linear_idx_start;
    unsigned int current_proj_split_size;
    size_t num_bytes_img_curr;
    size_t img_linear_idx_start;
    for( unsigned int proj=0;proj<split_projections;proj++){
        // Crop the last one, as its likely its not completely divisible.
        current_proj_split_size=(nalpha+split_projections-1)/split_projections;
        current_proj_split_size=((proj+1)*current_proj_split_size<nalpha)?  current_proj_split_size:  nalpha-current_proj_split_size*proj;
        
        //Get the linear index where the current memory chunk starts.
        proj_linear_idx_start=(unsigned long long)((nalpha+split_projections-1)/split_projections)*(unsigned long long)proj*(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV;
        
        // Now get the projections on memory
        cudaTextureObject_t *texProj = new cudaTextureObject_t[deviceCount];
        cudaArray **d_cuArrTex = new cudaArray*[deviceCount];
        CreateTexture2Spherical(deviceCount,&projections[proj_linear_idx_start],geo,d_cuArrTex,current_proj_split_size,texProj);
        
        
        for(unsigned int img_slice=0;img_slice<split_image;img_slice++){
            if(proj==0){
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    cudaMemset(dimage[dev],0,num_bytes_img);
                    cudaCheckErrors("memset fail");
                }
            }
            // If we have more than one chunck of projections, then we need to put back the previous results into the GPU
            // Exception: when each GPU is handling a single image chunk, we can leave it there.
            if(proj>0 && split_image>1) {
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    num_bytes_img_curr=geoArray[img_slice*deviceCount+dev].nVoxelX*geoArray[img_slice*deviceCount+dev].nVoxelY*geoArray[img_slice*deviceCount+dev].nVoxelZ*sizeof(float);
                    img_linear_idx_start=geo.nVoxelX*geo.nVoxelY*geoArray[0].nVoxelZ*(img_slice*deviceCount+dev);
                    cudaMemcpy(dimage[dev],&result[img_linear_idx_start], num_bytes_img_curr, cudaMemcpyHostToDevice);
                    cudaCheckErrors("cudaMemcpy previous result fail");
                }
                
            }
            
            for (dev = 0; dev < deviceCount; dev++){
                //Safety:
                // Depends on the amount of GPUs, the case where a image slice is zero hight can happen.
                // Just break the loop if we reached that point
                if(geoArray[img_slice*deviceCount+dev].nVoxelZ==0)
                    break;
                
                cudaSetDevice(dev);
                
                
                int divx,divy,divz;
        
                divx=32;
                divy=32;
                divz=1;
                dim3 grid((geo.nVoxelX+divx-1)/divx,
                        (geo.nVoxelY+divy-1)/divy,
                        (geoArray[img_slice*deviceCount+dev].nVoxelZ+divz-1)/divz);
                dim3 block(divx,divy,divz);
                Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, offDetec,source;

                // Main loop
                for (unsigned int i=0;i<current_proj_split_size;i++){
                    unsigned int currProjNumber_global=i+proj*(nalpha+split_projections-1)/split_projections;
                    geo.alpha=-angles[currProjNumber_global*3];
                    geo.theta=-angles[currProjNumber_global*3+1];
                    geo.psi  =-angles[currProjNumber_global*3+2];
                    
                    computeDeltasCubeSpherical(geoArray[img_slice*deviceCount+dev],currProjNumber_global,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
                    
                    offOrig.x=geo.offOrigX[currProjNumber_global];
                    offOrig.y=geo.offOrigY[currProjNumber_global];
                    offOrig.z=geoArray[img_slice*deviceCount+dev].offOrigZ[currProjNumber_global];
                    offDetec.x=geo.offDetecU[currProjNumber_global];
                    offDetec.y=geo.offDetecV[currProjNumber_global];
                    
                    kernelPixelBackprojectionSpherical<<<grid,block>>>
                            (geoArray[img_slice*deviceCount+dev],dimage[dev],currProjNumber_global,geo.COR[currProjNumber_global],geo.DSD[currProjNumber_global],geo.DSO[currProjNumber_global],cos(geo.alpha),sin(geo.alpha),deltaX,deltaY,deltaZ,xyzOrigin,offOrig,offDetec,source,texProj[dev]);
                    cudaCheckErrors("Kernel fail");
                }
                
            }
            
            //If this is the last time this image chiunk is going to be backprojected, then lets apply the weigth.
            if (proj==split_projections-1){
                cudaDeviceSynchronize();
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    matrixConstantMultiply_sp<<<60,MAXTREADS>>>( geoArray[img_slice*deviceCount+dev],dimage[dev],geo.dVoxelX*geo.dVoxelY*geo.dVoxelZ/(geo.dDetecU*geo.dDetecV));
                }
            }
            
            
            //////////////////////////////////////////////////////////////////////////////////////
            // END Main reconstruction loop: go through projections (rotation angles) and backproject
            //////////////////////////////////////////////////////////////////////////////////////
            // Now that we have backprojected all the image chunk (for the current projection chunk)
            // we need to take it out of the GPU
            // Exception: when each GPU is handlin a single image chunk and we have not finished with the projections
            cudaDeviceSynchronize();
            if(proj==split_projections-1 || split_image>1){
                
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    num_bytes_img_curr=geoArray[img_slice*deviceCount+dev].nVoxelX*geoArray[img_slice*deviceCount+dev].nVoxelY*geoArray[img_slice*deviceCount+dev].nVoxelZ*sizeof(float);
                    img_linear_idx_start=geo.nVoxelX*geo.nVoxelY*geoArray[0].nVoxelZ*(img_slice*deviceCount+dev);
                    cudaMemcpy(&result[img_linear_idx_start], dimage[dev], num_bytes_img_curr, cudaMemcpyDeviceToHost);
                    cudaCheckErrors("cudaMemcpy result fail");
                }
            }
            cudaDeviceSynchronize();
        }
        
        //Before the next iteration of projection slices, delete the current textures.
        
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(dev);
            cudaDestroyTextureObject(texProj[dev]);
            cudaFreeArray(d_cuArrTex[dev]);
            
        }
    }
    
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaFree(dimage[dev]);
    }
    freeGeoArray(split_image*deviceCount,geoArray);
    cudaCheckErrors("cudaFree d_imagedata fail");
    cudaDeviceReset(); // For the Nvidia Visual Profiler
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
void CreateTexture2Spherical(int num_devices,const float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage)
{
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    for (unsigned int i = 0; i < num_devices; i++){
        cudaSetDevice(i);
        
        //cudaArray Descriptor
        const cudaExtent extent = make_cudaExtent(geo.nDetecU, geo.nDetecV, nangles);
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        //cuda Array
        cudaMalloc3DArray(&d_cuArrTex[i], &channelDesc, extent);
        cudaCheckErrors("Texture memory allocation fail");
        cudaMemcpy3DParms copyParams = {0};
        
        
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)projectiondata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[i];
        copyParams.extent   = extent;
        copyParams.kind     = cudaMemcpyHostToDevice;
        cudaMemcpy3D(&copyParams);
        cudaCheckErrors("Texture memory data copy fail");
        //Array creation End
        
        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArrTex[i];
        cudaTextureDesc     texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = false;
        texDescr.filterMode = cudaFilterModeLinear;
        texDescr.addressMode[0] = cudaAddressModeBorder;
        texDescr.addressMode[1] = cudaAddressModeBorder;
        texDescr.addressMode[2] = cudaAddressModeBorder;
        texDescr.readMode = cudaReadModeElementType;
        cudaCreateTextureObject(&texImage[i], &texRes, &texDescr, NULL);
        cudaCheckErrors("Texture object creation fail");
    }
}
#ifndef BACKPROJECTION_HPP

void checkDevices(void){
    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning thrown)
    int dev;
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    char * devicenames;
    cudaDeviceProp deviceProp;
    for (dev = 0; dev < deviceCount; dev++) {
        cudaSetDevice(dev);
        cudaGetDeviceProperties(&deviceProp, dev);
        if (dev>0){
            if (strcmp(devicenames,deviceProp.name)!=0){
                mexWarnMsgIdAndTxt("Atb:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed. \n Siddon_projection.cu line 275.");
                break;
            }
        }
        devicenames=deviceProp.name;
    }
}
void splitCTbackprojection(int deviceCount,Geometry geo,int nalpha, unsigned int* split_image, unsigned int * split_projections){
    
    
    size_t memfree;
    size_t memtotal;
    size_t mem_GPU_global;
    for (unsigned int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(dev);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            mexErrMsgIdAndTxt("minimizeTV:POCS_TV:GPU","One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n Free the GPU to run TIGRE\n");
        }
        cudaCheckErrors("Check mem error");
        
        mem_GPU_global=(memfree<mem_GPU_global)?memfree:mem_GPU_global;
    }
    mem_GPU_global=(size_t)((double)mem_GPU_global*0.95); // lets leave 10% for the GPU. Too much? maybe, but probably worth saving.

    // Compute how much memory each of the relevant memory pieces need
    size_t mem_image=       (unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
    size_t mem_image_slice= (unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)VOXELS_PER_THREAD*sizeof(float);
    size_t mem_proj=        (unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV*sizeof(float);
    
    // Initialize variables for spliting procedure choosing algorithm.
    
    
    // Does everything fit in the GPU?
    
    if(mem_image+mem_proj*nalpha<mem_GPU_global){
        // We only need to split if we have extra GPUs
        *split_image=1;
        *split_projections=1;
    }
    // We know we need to split, but:
    // Do all projections fit on the GPU (with some slack for the image)??
    else if(mem_proj*nalpha+mem_image_slice <mem_GPU_global){
        // We should then store all the projections in GPU and split the image backprojection in as big chunks as possible
        *split_projections=1;
        size_t mem_free=mem_GPU_global-nalpha*mem_proj;
        // How many slices can we fit on the free memory we have? We'd like to keep these slices full as it increases
        // the kernels performance to run image chuncks that are VOXELS_PER_THREAD in z.
        unsigned int total_slices_img=(geo.nVoxelZ+VOXELS_PER_THREAD-1)/VOXELS_PER_THREAD;
        // Split:
        // mem_free/mem_image_slice == how many slices fit in each GPU
        // total_slices_img/deviceCount == How many slices we need each GPU to evaluate.
        *split_image=(unsigned int)ceil((float)total_slices_img/(float)deviceCount/(float)(mem_free/mem_image_slice));
        
    }
    // They do not fit in memory. We need to split both projections and images. Its OK, we'll survive.
    else
    {
        // Because I think memory transfer is the slowest operation if broken apart (TODO test this assumtion) we should minimize the amount
        // of those that happen. So the target split is that one that requires less copys of projection data to memory.
        // Lets assume to start with that we only need 1 slice of image in memory, the rest is fo rprojections
        size_t mem_free=mem_GPU_global-mem_image_slice;
        *split_projections=(mem_proj*nalpha+mem_free-1)/mem_free;
        // Now knowing how many splits we have for projections, we can recompute how many slices of image actually
        // fit on the GPU. Must be more than 0 obviously.
        mem_free=mem_GPU_global-mem_proj*nalpha/(*split_projections);
        unsigned int total_slices_img=(geo.nVoxelZ+VOXELS_PER_THREAD-1)/VOXELS_PER_THREAD;
        // Split:
        // mem_free/mem_image_slice == how many slices fit in each GPU
        // total_slices_img/deviceCount == How many slices we need each GPU to evaluate.
        *split_image=(unsigned int)ceil((float)total_slices_img/(float)deviceCount/(float)(mem_free/mem_image_slice));
        
    }
}
//______________________________________________________________________________
//
//      Function:       createGeoArray
//
//      Description:    This code generates the geometries needed to split the image properly in
//                      cases where the entire image does not fit in the memory of the GPU
//______________________________________________________________________________

void createGeoArray(unsigned int image_splits, Geometry geo,Geometry* geoArray, unsigned int nangles){
    
    
    unsigned int  splitsize=(geo.nVoxelZ+image_splits-1)/image_splits;
    
    for(unsigned int sp=0;sp<image_splits;sp++){
        geoArray[sp]=geo;
        // All of them are splitsize, but the last one, possible
        geoArray[sp].nVoxelZ=((sp+1)*splitsize<geo.nVoxelZ)?  splitsize:  max(geo.nVoxelZ-splitsize*sp,0);
        geoArray[sp].sVoxelZ= geoArray[sp].nVoxelZ* geoArray[sp].dVoxelZ;
        
        // We need to redefine the offsets, as now each subimage is not aligned in the origin.
        geoArray[sp].offOrigZ=(float *)malloc(nangles*sizeof(float));
        for (unsigned int i=0;i<nangles;i++){
            geoArray[sp].offOrigZ[i]=geo.offOrigZ[i]-geo.sVoxelZ/2+sp*geoArray[0].sVoxelZ+geoArray[sp].sVoxelZ/2;
        }
    }
    
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
#endif