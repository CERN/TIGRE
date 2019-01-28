/*-------------------------------------------------------------------------
 *
 * CUDA function for backrpojection using FDK weigts for CBCT
 *
 *
 * CODE by  Ander Biguri
 *          Optimized and modified by RB
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

#define  PI_2 1.57079632679489661923
#include <algorithm>
#include <cuda_runtime_api.h>
#include <cuda.h>
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
    
// this definitionmust go here.
    void CreateTexture2(int num_devices,const float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage);



__global__ void matrixConstantMultiply(const Geometry geo,float* image,float constant){
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    for(; idx<geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ; idx+=gridDim.x*blockDim.x) {
        image[idx]*=constant;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RB, 10/31/2016: Add constant memory arrays to store parameters for all projections to be analyzed during a single kernel call
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The optimal values of two constants obtained by RB on NVIDIA Quadro K2200 (4 GB RAM, 640 CUDA cores) for 512^3 volume and 512^3 projections (512 proj, each 512 x 512) were:
// PROJ_PER_KERNEL = 32 or 16 (very similar times)
// VOXELS_PER_THREAD = 8
// Speedup of the entire FDK backprojection (not only kernel run, also memcpy etc.) was nearly 4x relative to the original (single projection, single voxel per thread) code.
// (e.g. 16.2 s vs. ~62 s).

const int PROJ_PER_KERNEL = 32;  // Number of 2D projections to be analyzed by a single thread. This can be tweaked to see what works best. 32 was the optimal value in the paper by Zinsser and Keck.
const int VOXELS_PER_THREAD = 8;  // Number of voxels to be computed by s single thread. Can be tweaked to see what works best. 4 was the optimal value in the paper by Zinsser and Keck.

// We have PROJ_PER_KERNEL projections and we need 6 parameters for each projection:
//   deltaX, deltaY, deltaZ, xyzOrigin, offOrig, offDetec
// So we need to keep PROJ_PER_KERNEL*6 values in our deltas array FOR EACH CALL to our main kernel
// (they will be updated in the main loop before each kernel call).

__constant__ Point3D projParamsArray2Dev[7*PROJ_PER_KERNEL];  // Dev means it is on device

// We also need a corresponding array on the host side to be filled before each kernel call, then copied to the device (array in constant memory above)
Point3D projParamsArray2Host[7*PROJ_PER_KERNEL];   // Host means it is host memory

// Now we also need to store sinAlpha and cosAlpha for each projection (two floats per projection)
__constant__ float projSinCosArray2Dev[5*PROJ_PER_KERNEL];

float projSinCosArray2Host[5*PROJ_PER_KERNEL];
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END RB, 10/31/2016: Add constant memory arrays to store parameters for all projections to be analyzed during a single kernel call
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
//
//      Function:       kernelPixelBackprojectionFDK
//
//      Description:    Main FDK backprojection kernel
//______________________________________________________________________________

__global__ void kernelPixelBackprojection(const Geometry geo, float* image,const int currProjSetNumber, const int totalNoOfProjections, cudaTextureObject_t tex)
{
    
    unsigned long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long indX = blockIdx.x * blockDim.x + threadIdx.x;
    // unsigned long startIndZ = blockIdx.z * blockDim.z + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    unsigned long startIndZ = blockIdx.z * VOXELS_PER_THREAD + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    //Make sure we dont go out of bounds
    if (indX>=geo.nVoxelX | indY>=geo.nVoxelY |startIndZ>=geo.nVoxelZ)
        return;
    
    // We'll keep a local auxiliary array of values of a column of voxels that this thread will update
    float voxelColumn[VOXELS_PER_THREAD];
    
    // First we need to copy the curent 3D volume values from the column to our auxiliary array so that we can then
    // work on them (update them by computing values from multiple projections) locally - avoiding main memory reads/writes
    
    int colIdx;
#pragma unroll
    for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
    {
        unsigned long indZ = startIndZ + colIdx;
        // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
        // be trying to copy the out of bounds values back to the 3D volume anyway (bounds checks will be done in the final loop where the updated values go back to the main volume)
        if(indZ>=geo.nVoxelZ)
            break;   // break the loop.
        
        unsigned long long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
        voxelColumn[colIdx] = image[idx];   // Read the current volume value that we'll update by computing values from MULTIPLE projections (not just one)
        // We'll be updating the local (register) variable, avoiding reads/writes from the slow main memory.
    }  // END copy 3D volume voxels to local array
    
    // Now iterate through projections
#pragma unroll
    for(int projNumber=0; projNumber<PROJ_PER_KERNEL; projNumber++)
    {
        // Get the current parameters from parameter arrays in constant memory.
        int indAlpha = currProjSetNumber*PROJ_PER_KERNEL+projNumber;  // This is the ABSOLUTE projection number in the projection array
        
        // Our currImageVal will be updated by hovewer many projections we had left in the "remainder" - that's OK.
        if(indAlpha>=totalNoOfProjections)
            break;
        
        Point3D deltaX = projParamsArray2Dev[7*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
        Point3D deltaY = projParamsArray2Dev[7*projNumber+1];
        Point3D deltaZ = projParamsArray2Dev[7*projNumber+2];
        Point3D xyzOrigin = projParamsArray2Dev[7*projNumber+3];
        Point3D xyzOffset = projParamsArray2Dev[7*projNumber+4];
        Point3D uv0Offset = projParamsArray2Dev[7*projNumber+5];
        Point3D S = projParamsArray2Dev[7*projNumber+6];
        
        float sinalpha = projSinCosArray2Dev[5*projNumber];     // 2*projNumber because we have 2 float (sin or cos angle) values per projection
        float cosalpha = projSinCosArray2Dev[5*projNumber+1];
        float COR = projSinCosArray2Dev[5*projNumber+2];
        float DSD = projSinCosArray2Dev[5*projNumber+3];
        float DSO = projSinCosArray2Dev[5*projNumber+4];
        // Precomputations for the weights:
        //Real coords of Source
        // We already have S.x (geo.DSO), and S.y and S.z are always zero. we just need to rotate
        Point3D realS;
        realS.x= DSO*cosalpha;
        realS.y=-DSO*sinalpha;
        realS.z=0;
        
        
        Point3D realvoxel_init;
        realvoxel_init.x=-geo.sVoxelX/2+geo.dVoxelX/2+xyzOffset.x;
        realvoxel_init.y=-geo.sVoxelY/2+geo.dVoxelY/2+xyzOffset.y;
        realvoxel_init.z=-geo.sVoxelZ/2+geo.dVoxelZ/2+xyzOffset.z;
        
        // Now iterate through Z in our voxel column FOR A GIVEN PROJECTION
#pragma unroll
        for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
        {
            unsigned long indZ = startIndZ + colIdx;
            
            // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
            // be trying to copy the out of bounds values anyway (bounds checks will be done in the final loop where the values go to the main volume)
            if(indZ>=geo.nVoxelZ)
                break;   // break the loop.
            
            // "XYZ" in the scaled coordinate system of the current point. The image is rotated with the projection angles.
            Point3D P;
            P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
            P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y)-COR/geo.dDetecU;
            P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
            
            // This is the vector defining the line from the source to the Voxel
            float vectX,vectY,vectZ;
            vectX=(P.x -S.x);
            vectY=(P.y -S.y);
            vectZ=(P.z -S.z);
            
            // Get the coordinates in the detector UV where the mid point of the voxel is projected.
            float t=(DSO-DSD /*-DOD*/ - S.x)/vectX;
            float y,z;
            y=vectY*t+S.y;
            z=vectZ*t+S.z;
            float u,v;
            u=y+geo.nDetecU/2;
            v=z+geo.nDetecV/2;
            float sample=tex3D<float>(tex, v, u ,indAlpha+0.5f);
            float weigth=0;
            //
            //
            //
            // IMPORTANT: The weights are almost 50% of the computational time. Is there a way of speeding this up??
            //
            //Real coordinates of Voxel. Instead of reverting the tranformation, its less math (faster) to compute it from the indexes.
            Point3D realvoxel;
            
            realvoxel.x=realvoxel_init.x+indX*geo.dVoxelX; 
            realvoxel.y=realvoxel_init.y+indY*geo.dVoxelY; 
            realvoxel.z=realvoxel_init.z+indZ*geo.dVoxelZ; 

            
            
            // Real XYZ coordinates of Detector.
            Point3D realD, realDaux;
            // We know the index of the detector (u,v). Start from there.
            realDaux.x=-(DSD-DSO);
            realDaux.y=-geo.sDetecU/2+geo.dDetecU/2 + u*geo.dDetecU +uv0Offset.x;
            realD.z   =-geo.sDetecV/2+geo.dDetecV/2 + v*geo.dDetecV +uv0Offset.y;
            //rotate the detector
            realD.x= realDaux.x*cosalpha  + realDaux.y*sinalpha; //sin(-x)=-sin(x) , cos(-x)=cos(x)
            realD.y=-realDaux.x*sinalpha  + realDaux.y*cosalpha; //sin(-x)=-sin(x) , cos(-x)=cos(x)
            float L=0,l=0;
            
            L = sqrt( (realS.x-realD.x)*(realS.x-realD.x)+ (realS.y-realD.y)*(realS.y-realD.y)+ (realD.z)*(realD.z)); // Sz=0 always.
            l = sqrt( (realS.x-realvoxel.x)*(realS.x-realvoxel.x)
            + (realS.y-realvoxel.y)*(realS.y-realvoxel.y)
            + (realS.z-realvoxel.z)*(realS.z-realvoxel.z));
            weigth=L*L*L/(DSD*l*l);
            
            // Get Value in the computed (U,V) and multiply by the corresponding weigth.
            // indAlpha is the ABSOLUTE number of projection in the projection array (NOT the current number of projection set!)
            voxelColumn[colIdx]+=sample* weigth;
        }  // END iterating through column of voxels
        
    }  // END iterating through multiple projections
    
    // And finally copy the updated local voxelColumn array back to our 3D volume (main memory)
#pragma unroll
    for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
    {
        unsigned long indZ = startIndZ + colIdx;
        // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
        // be trying to copy the out of bounds values back to the 3D volume anyway (bounds checks will be done in the final loop where the values go to the main volume)
        if(indZ>=geo.nVoxelZ)
            break;   // break the loop.
        
        unsigned long long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
        image[idx] = voxelColumn[colIdx];   // Read the current volume value that we'll update by computing values from MULTIPLE projections (not just one)
        // We'll be updating the local (register) variable, avoiding reads/writes from the slow main memory.
        // According to references (Papenhausen), doing = is better than +=, since += requires main memory read followed by a write.
        // We did all the reads into the local array at the BEGINNING of this kernel. According to Papenhausen, this type of read-write split is
        // better for avoiding memory congestion.
    }  // END copy updated voxels from local array to our 3D volume
    
}  // END kernelPixelBackprojectionFDK




//______________________________________________________________________________
//
//      Function:       voxel_backprojection
//
//      Description:    Main host function for FDK backprojection (invokes the kernel)
//______________________________________________________________________________

int voxel_backprojection2(float const * const projections, Geometry geo, float* result,float const * const alphas, int nalpha){
    
    
    
    
    // Prepare for MultiGPU
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    cudaCheckErrors("Device query fail");
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Atb:Voxel_backprojection:GPUselect","There are no available device(s) that support CUDA\n");
    }
    
    int dev;
    checkDevices();
    
    
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
        CreateTexture2(deviceCount,&projections[proj_linear_idx_start],geo,d_cuArrTex,current_proj_split_size,texProj);
        
        
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
                
                // RB: Use the optimal (in their tests) block size from paper by Zinsser and Keck (16 in x and 32 in y).
                // I tried different sizes and shapes of blocks (tiles), but it does not appear to significantly affect trhoughput, so
                // let's stick with the values from Zinsser and Keck.
                divx=16;
                divy=32;
                divz=VOXELS_PER_THREAD;      // We now only have 32 x 16 threads per block (flat tile, see below), BUT each thread works on a Z column of VOXELS_PER_THREAD voxels, so we effectively need fewer blocks!
                dim3 grid((geo.nVoxelX+divx-1)/divx,
                        (geo.nVoxelY+divy-1)/divy,
                        (geoArray[img_slice*deviceCount+dev].nVoxelZ+divz-1)/divz);
                
                dim3 block(divx,divy,1);    // Note that we have 1 in the Z size, not divz, since each thread works on a vertical set of VOXELS_PER_THREAD voxels (so we only need a "flat" tile of threads, with depth of 1)
                
                //////////////////////////////////////////////////////////////////////////////////////
                // Main reconstruction loop: go through projections (rotation angles) and backproject
                //////////////////////////////////////////////////////////////////////////////////////
                
                // Since we'll have multiple projections processed by a SINGLE kernel call, compute how many
                // kernel calls we'll need altogether.
                
                
                int noOfKernelCalls = (current_proj_split_size+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_KERNEL
                
                for (unsigned int i=0; i<noOfKernelCalls; i++){
                    // Now we need to generate and copy all data for PROJ_PER_KERNEL projections to constant memory so that our kernel can use it
                    int j;
                    for(j=0; j<PROJ_PER_KERNEL; j++){
                        unsigned int currProjNumber_slice=i*PROJ_PER_KERNEL+j;
                        unsigned int currProjNumber_global=i*PROJ_PER_KERNEL+j+proj*(nalpha+split_projections-1)/split_projections;
                        
                        if(currProjNumber_slice>=current_proj_split_size)
                            break;  // Exit the loop. Even when we leave the param arrays only partially filled, this is OK, since the kernel will check bounds anyway.
                        
                        Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, offDetec,source;
                        float sinalpha,cosalpha;
                        
                        geo.alpha=-alphas[currProjNumber_global*3];
                        sinalpha=sin(geo.alpha);
                        cosalpha=cos(geo.alpha);
                        
                        projSinCosArray2Host[5*j]=sinalpha;  // 2*j because we have 2 float (sin or cos angle) values per projection
                        projSinCosArray2Host[5*j+1]=cosalpha;
                        projSinCosArray2Host[5*j+2]=geo.COR[currProjNumber_global];
                        projSinCosArray2Host[5*j+3]=geo.DSD[currProjNumber_global];
                        projSinCosArray2Host[5*j+4]=geo.DSO[currProjNumber_global];
                        
                        computeDeltasCube(geoArray[img_slice*deviceCount+dev],geo.alpha,currProjNumber_global,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
                        
                        offOrig.x=geo.offOrigX[currProjNumber_global];
                        offOrig.y=geo.offOrigY[currProjNumber_global];
                        offOrig.z=geoArray[img_slice*deviceCount+dev].offOrigZ[currProjNumber_global];
                        
                        offDetec.x=geo.offDetecU[currProjNumber_global];
                        offDetec.y=geo.offDetecV[currProjNumber_global];
                        
                        projParamsArray2Host[7*j]  =deltaX;		// 7*j because we have 7 Point3D values per projection
                        projParamsArray2Host[7*j+1]=deltaY;
                        projParamsArray2Host[7*j+2]=deltaZ;
                        projParamsArray2Host[7*j+3]=xyzOrigin;
                        projParamsArray2Host[7*j+4]=offOrig;
                        projParamsArray2Host[7*j+5]=offDetec;
                        projParamsArray2Host[7*j+6]=source;
                        
                    }   // END for (preparing params for kernel call)
                    
                    // Copy the prepared parameter arrays to constant memory to make it available for the kernel
                    cudaMemcpyToSymbolAsync(projSinCosArray2Dev, projSinCosArray2Host, sizeof(float)*5*PROJ_PER_KERNEL);
                    cudaMemcpyToSymbolAsync(projParamsArray2Dev, projParamsArray2Host, sizeof(Point3D)*7*PROJ_PER_KERNEL);
                    
                    kernelPixelBackprojection<<<grid,block>>>(geoArray[img_slice*deviceCount+dev],dimage[dev],i,current_proj_split_size,texProj[dev]);
                    
                    cudaCheckErrors("Kernel fail");
                }  // END for
            }
            
            //If this is the last time this image chiunk is going to be backprojected, then lets apply the weigth.
            if (proj==split_projections-1){
                cudaDeviceSynchronize();
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(dev);
                    matrixConstantMultiply<<<60,MAXTREADS>>>( geoArray[img_slice*deviceCount+dev],dimage[dev],geo.dVoxelX*geo.dVoxelY*geo.dVoxelZ/(geo.dDetecU*geo.dDetecV));
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
    
}  // END voxel_backprojection



void CreateTexture2(int num_devices,const float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage)
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
    mem_GPU_global=(size_t)((double)mem_GPU_global*0.95); // lets leave 5% for the GPU. Too much? maybe, but probably worth saving.
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
void computeDeltasCube(Geometry geo, float alpha,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S)
{
    Point3D P0, Px0,Py0,Pz0, source;
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
    
    //detector offset
    P.z =P.z-geo.offDetecV[i];          P.y =P.y-geo.offDetecU[i];
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
    
    source.x=geo.DSD[i]; //allready offseted for rotation of teh detector
    source.y=-geo.offDetecU[i];
    source.z=-geo.offDetecV[i];
    rollPitchYawT(geo,i,&source);
    
    source.x=source.x-(geo.DSD[i]-geo.DSO[i]);
    
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
}  // END computeDeltasCube

void rollPitchYawT(Geometry geo,int i, Point3D* point){
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x=cos(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.x
            +sin(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.y
            -sin(geo.dPitch[i])*auxPoint.z;
    
    
    point->y=(cos(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) - sin(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.x
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) + cos(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.y
            +cos(geo.dPitch[i])*sin(geo.dYaw[i])*auxPoint.z;
    
    
    point->z=(cos(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) + sin(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.x
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) - cos(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.y
            +cos(geo.dPitch[i])*cos(geo.dYaw[i])*auxPoint.z;
    
}
#endif