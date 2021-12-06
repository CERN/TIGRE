/*-------------------------------------------------------------------------
 *
 * CUDA function for backrpojection  for parallel beam
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
#include "voxel_backprojection.hpp"
#include "voxel_backprojection_parallel.hpp"

#include "TIGRE_common.hpp"
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
void CreateTextureParallel( float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream, bool allocate);

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

__constant__ Point3D projParamsArrayDevParallel[6*PROJ_PER_KERNEL];  // Dev means it is on device

// We also need a corresponding array on the host side to be filled before each kernel call, then copied to the device (array in constant memory above)
// Point3D projParamsArrayHostParallel[6*PROJ_PER_KERNEL];   // Host means it is host memory

// Now we also need to store sinAlpha and cosAlpha for each projection (two floats per projection)
__constant__ float projSinCosArrayDevParallel[3*PROJ_PER_KERNEL];

// float projSinCosArrayHostParallel[3*PROJ_PER_KERNEL];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END RB, 10/31/2016: Add constant memory arrays to store parameters for all projections to be analyzed during a single kernel call
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//______________________________________________________________________________
//
//      Function:       kernelPixelBackprojectionFDK
//
//      Description:    Main FDK backprojection kernel
//______________________________________________________________________________

__global__ void kernelPixelBackprojection_parallel(const Geometry geo, float* image,const int currProjSetNumber, const int totalNoOfProjections,cudaTextureObject_t tex)
{
    
    // Old kernel call signature:
    // kernelPixelBackprojectionFDK<<<grid,block>>>(geo,dimage,i,deltaX,deltaY,deltaZ,xyzOrigin,offOrig,offDetec,sinalpha,cosalpha);
    // We just read in most of the params from the constant memory instead of getting them from the param list.
    // This is because we now have MANY params, since single kernel processes more than one projection!
    /* __global__ void kernelPixelBackprojectionFDK(const Geometry geo,
     * float* image,
     * const int indAlpha,
     * const Point3D deltaX ,
     * const Point3D deltaY,
     * const Point3D deltaZ,
     * const Point3D xyzOrigin,
     * const Point3D xyzOffset,
     * const Point3D uv0Offset,
     * const float sinalpha,
     * const float cosalpha){
     */
    unsigned long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long indX = blockIdx.x * blockDim.x + threadIdx.x;
    // unsigned long startIndZ = blockIdx.z * blockDim.z + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    unsigned long startIndZ = blockIdx.z * VOXELS_PER_THREAD + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    //Make sure we don't go out of bounds
    if (indX>=geo.nVoxelX || indY>=geo.nVoxelY || startIndZ>=geo.nVoxelZ)
        return;
    
    // We'll keep a local auxiliary array of values of a column of voxels that this thread will update
    float voxelColumn[VOXELS_PER_THREAD];
    
    // First we need to copy the curent 3D volume values from the column to our auxiliary array so that we can then
    // work on them (update them by computing values from multiple projections) locally - avoiding main memory reads/writes
    
    int colIdx;
    
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
    for(int projNumber=0; projNumber<PROJ_PER_KERNEL; projNumber++)
    {
        // Get the current parameters from parameter arrays in constant memory.
        int indAlpha = currProjSetNumber*PROJ_PER_KERNEL+projNumber;  // This is the ABSOLUTE projection number in the projection array
        
        // Our currImageVal will be updated by hovewer many projections we had left in the "remainder" - that's OK.
        if(indAlpha>=totalNoOfProjections)
            break;
        
        Point3D deltaX = projParamsArrayDevParallel[6*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
        Point3D deltaY = projParamsArrayDevParallel[6*projNumber+1];
        Point3D deltaZ = projParamsArrayDevParallel[6*projNumber+2];
        Point3D xyzOrigin = projParamsArrayDevParallel[6*projNumber+3];
        Point3D xyzOffset = projParamsArrayDevParallel[6*projNumber+4];
        Point3D S = projParamsArrayDevParallel[6*projNumber+5];
        
        float DSD = projSinCosArrayDevParallel[3*projNumber];     // 2*projNumber because we have 2 float (sin or cos angle) values per projection
        float DSO = projSinCosArrayDevParallel[3*projNumber+1];
        float COR = projSinCosArrayDevParallel[3*projNumber+2];
        
        // Geometric trasnformations:
        //Source, scaled XYZ coordinates
        
        // Now iterate through Z in our voxel column FOR A GIVEN PROJECTION
        for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
        {
            unsigned long indZ = startIndZ + colIdx;
            
            // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
            // be trying to copy the out of bounds values anyway (bounds checks will be done in the final loop where the values go to the main volume)
            if(indZ>=geo.nVoxelZ)
                break;   // break the loop.
            
            // "XYZ" in the scaled coordinate system of the current point. The image is rotated with the projection angles.
            Point3D P;
            S.x=DSO;
            P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
            P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y)-COR/geo.dDetecU;
            P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
            S.y=P.y;S.z=P.z;
            
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
            u=y+geo.nDetecU/2.0f-0.5f;
            v=z+geo.nDetecV/2.0f-0.5f;
            
            
            
            // Get Value in the computed (U,V) and multiply by the corresponding weigth.
            // indAlpha is the ABSOLUTE number of projection in the projection array (NOT the current number of projection set!)
#if IS_FOR_MATLAB_TIGRE
            voxelColumn[colIdx]+=tex3D<float>(tex, v+0.5f, u+0.5f ,indAlpha+0.5f);
#else
            voxelColumn[colIdx]+=tex3D<float>(tex, u+0.5f, v+0.5f ,indAlpha+0.5f);
#endif
            
        }  // END iterating through column of voxels
        
    }  // END iterating through multiple projections
    
    // And finally copy the updated local voxelColumn array back to our 3D volume (main memory)
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
//      Function:       voxel_backprojection_parallel
//
//      Description:    Main host function for FDK backprojection (invokes the kernel)
//______________________________________________________________________________

int voxel_backprojection_parallel(float  *  projections, Geometry geo, float* result,float const * const alphas, int nalpha, const GpuIds& gpuids)
{
    if (gpuids.GetLength() == 0) {
        cudaSetDevice(0);
    } else {
        cudaSetDevice(gpuids[0]);
    }
    
    /*
     * Allocate texture memory on the device
     */
    // copy data to CUDA memory
    //If it is the first time, lets make sure our image is zeroed.
    int nStreamDevice=2;
    int nStreams=nStreamDevice;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));;
    
    for (int i = 0; i < nStreamDevice; ++i){
        cudaStreamCreate(&stream[i]);
        
        
    }
    //Pagelock memory for synchronous copy.
    // Lets try to make the host memory pinned:
    // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
    if (isHostRegisterSupported){
        cudaHostRegister(projections, (size_t)geo.nDetecU*(size_t)geo.nDetecV*(size_t)nalpha*(size_t)sizeof(float),cudaHostRegisterPortable);
    }
    cudaCheckErrors("Error pinning memory");
    
    
    // Allocate result image memory
    size_t num_bytes = geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ * sizeof(float);
    float* dimage;
    cudaMalloc((void**)&dimage, num_bytes);
    cudaMemset(dimage,0,num_bytes);
    cudaCheckErrors("cudaMalloc fail");
    
    
    Point3D* projParamsArrayHostParallel;
    cudaMallocHost((void**)&projParamsArrayHostParallel,6*PROJ_PER_KERNEL*sizeof(Point3D));
    float* projSinCosArrayHostParallel;
    cudaMallocHost((void**)&projSinCosArrayHostParallel,3*PROJ_PER_KERNEL*sizeof(float));
    
    
    // Texture buffer objects
    cudaTextureObject_t *texProj;
    cudaArray **d_cuArrTex;
    texProj =(cudaTextureObject_t*)malloc(2*sizeof(cudaTextureObject_t));
    d_cuArrTex =(cudaArray**)malloc(2*sizeof(cudaArray*));

    
    
    unsigned int proj_split_overlap_number;
    unsigned int split_projections=1;
    // Start with the main loop. The Projection data needs to be allocated and dealocated in the main loop
    // as due to the nature of cudaArrays, we can not reuse them. This should not be a problem for the fast execution
    // of the code, as repeated allocation and deallocation only happens when the projection data is very very big,
    // and therefore allcoation time should be negligible, fluctuation of other computations should mask the time.
    unsigned long long proj_linear_idx_start;
    unsigned int current_proj_split_size,current_proj_overlap_split_size;
    size_t num_bytes_img_curr;
    size_t img_linear_idx_start;
    
    
    current_proj_split_size=nalpha;
    // We are going to split it in the same amount of kernels we need to execute.
    proj_split_overlap_number=(current_proj_split_size+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;
    
    
    // Create pointer to pointers of projections and precompute their location and size.
    
    float ** partial_projection=(float**)malloc(current_proj_split_size*sizeof(float*));
    size_t * proj_split_size=(size_t*)malloc(current_proj_split_size*sizeof(size_t*));
    
    for(unsigned int proj_block_split=0; proj_block_split<proj_split_overlap_number;proj_block_split++){
        // Crop the last one, as its likely its not completely divisible.
        // now lets split this for simultanoeus memcopy and compute.
        // We want to make sure that if we can, we run PROJ_PER_KERNEL projections, to maximize kernel acceleration
        // current_proj_overlap_split_size units = angles
        current_proj_overlap_split_size=max((current_proj_split_size+proj_split_overlap_number-1)/proj_split_overlap_number,PROJ_PER_KERNEL);
        current_proj_overlap_split_size=(proj_block_split<proj_split_overlap_number-1)?current_proj_overlap_split_size:current_proj_split_size-(proj_split_overlap_number-1)*current_proj_overlap_split_size;
        //Get the linear index where the current memory chunk starts.
        
        proj_linear_idx_start=proj_block_split*max((current_proj_split_size+proj_split_overlap_number-1)/proj_split_overlap_number,PROJ_PER_KERNEL)*(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV;
        //Store result
        proj_split_size[proj_block_split]=current_proj_overlap_split_size;
        partial_projection[proj_block_split]=&projections[proj_linear_idx_start];
        
    }
    for(unsigned int proj_block_split=0; proj_block_split<proj_split_overlap_number;proj_block_split++){
        
        // Now get the projections on memory
        
        CreateTextureParallel(partial_projection[proj_block_split],geo,
                &d_cuArrTex[(proj_block_split%2)],
                proj_split_size[proj_block_split],
                &texProj   [(proj_block_split%2)],
                stream,
                (proj_block_split<2));// Only allocate if its the first 2 calls
        
  
        cudaStreamSynchronize(stream[0+1]);
        
        

        int divx,divy,divz;
        
        // RB: Use the optimal (in their tests) block size from paper by Zinsser and Keck (16 in x and 32 in y).
        // I tried different sizes and shapes of blocks (tiles), but it does not appear to significantly affect trhoughput, so
        // let's stick with the values from Zinsser and Keck.
        divx=16;
        divy=32;
        divz=VOXELS_PER_THREAD;      // We now only have 32 x 16 threads per block (flat tile, see below), BUT each thread works on a Z column of VOXELS_PER_THREAD voxels, so we effectively need fewer blocks!
        dim3 grid((geo.nVoxelX+divx-1)/divx,
                (geo.nVoxelY+divy-1)/divy,
                (geo.nVoxelZ+divz-1)/divz);
        
        dim3 block(divx,divy,1);    // Note that we have 1 in the Z size, not divz, since each thread works on a vertical set of VOXELS_PER_THREAD voxels (so we only need a "flat" tile of threads, with depth of 1)
        
        
        
        
        //////////////////////////////////////////////////////////////////////////////////////
        // Main reconstruction loop: go through projections (rotation angles) and backproject
        //////////////////////////////////////////////////////////////////////////////////////
        
        // Since we'll have multiple projections processed by a SINGLE kernel call, compute how many
        // kernel calls we'll need altogether.
        int noOfKernelCalls = (proj_split_size[proj_block_split]+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_KERNEL
        for (unsigned int i=0; i<noOfKernelCalls; i++)
        {
            // Now we need to generate and copy all data for PROJ_PER_KERNEL projections to constant memory so that our kernel can use it
            int j;
            for(j=0; j<PROJ_PER_KERNEL; j++)
            {
                int currProjNumber=i*PROJ_PER_KERNEL+j;
                unsigned int currProjNumber_slice=i*PROJ_PER_KERNEL+j;
                unsigned int currProjNumber_global=i*PROJ_PER_KERNEL+j                                                                          // index within kernel
                        +proj_block_split*max(current_proj_split_size/proj_split_overlap_number,PROJ_PER_KERNEL); // indexof overlap current split
                if(currProjNumber_slice>=proj_split_size[proj_block_split])
                    break;  // Exit the loop. Even when we leave the param arrays only partially filled, this is OK, since the kernel will check bounds anyway.
                
                if(currProjNumber_global>=nalpha)
                    break;  // Exit the loop. Even when we leave the param arrays only partially filled, this is OK, since the kernel will check bounds anyway.
                
                Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, /*offDetec,*/source;
                float sinalpha,cosalpha;
                
                geo.alpha=-alphas[currProjNumber_global*3];
                geo.theta=-alphas[currProjNumber_global*3+1];
                geo.psi  =-alphas[currProjNumber_global*3+2];
                
                //sinalpha=sin(geo.alpha);
//            cosalpha=cos(geo.alpha);
                
                projSinCosArrayHostParallel[3*j]=geo.DSD[currProjNumber_global];  // 3*j because we have 3 float (sin or cos angle) values per projection
                projSinCosArrayHostParallel[3*j+1]=geo.DSO[currProjNumber_global];
                projSinCosArrayHostParallel[3*j+2]=geo.COR[currProjNumber_global];
                
                //computeDeltasCubeParallel(geo,geo.alpha,currProjNumber,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
                computeDeltasCubeParallel(geo,currProjNumber_global,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
                
                offOrig.x=geo.offOrigX[currProjNumber_global];
                offOrig.y=geo.offOrigY[currProjNumber_global];
                
                
                projParamsArrayHostParallel[6*j]=deltaX;		// 6*j because we have 6 Point3D values per projection
                projParamsArrayHostParallel[6*j+1]=deltaY;
                projParamsArrayHostParallel[6*j+2]=deltaZ;
                projParamsArrayHostParallel[6*j+3]=xyzOrigin;
                projParamsArrayHostParallel[6*j+4]=offOrig;
                projParamsArrayHostParallel[6*j+5]=source;
            }   // END for (preparing params for kernel call)
            
            // Copy the prepared parameter arrays to constant memory to make it available for the kernel
            
            cudaMemcpyToSymbolAsync(projSinCosArrayDevParallel, projSinCosArrayHostParallel, sizeof(float)*3*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[0]);
            cudaMemcpyToSymbolAsync(projParamsArrayDevParallel, projParamsArrayHostParallel, sizeof(Point3D)*6*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[0]);
            cudaStreamSynchronize(stream[0]);

            kernelPixelBackprojection_parallel<<<grid,block,0,stream[0]>>>(geo,dimage,i,proj_split_size[proj_block_split],texProj[(proj_block_split%2)]);
        }  // END for
        
        //////////////////////////////////////////////////////////////////////////////////////
        // END Main reconstruction loop: go through projections (rotation angles) and backproject
        //////////////////////////////////////////////////////////////////////////////////////
    }
    cudaDeviceSynchronize();
    cudaMemcpy(result, dimage, num_bytes, cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy result fail");
    
    free(partial_projection);
    free(proj_split_size);
        
  bool two_buffers_used=((((nalpha+split_projections-1)/split_projections)+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL)>1;
  for(unsigned int i=0; i<2;i++){ // 2 buffers (if needed, maybe only 1)
      if (!two_buffers_used && i==1)
          break;            
          cudaDestroyTextureObject(texProj[i]);
          cudaFreeArray(d_cuArrTex[i]);
  }
    cudaFreeHost(projSinCosArrayHostParallel);
    cudaFreeHost(projParamsArrayHostParallel);
    
    cudaFree(dimage);
     if (isHostRegisterSupported){
        cudaHostUnregister(projections);
    }
    for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]);

//     cudaDeviceReset();
    return 0;
    
}  // END voxel_backprojection

void computeDeltasCubeParallel(Geometry geo, int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D *S)
{
    
    Point3D P, Px,Py,Pz;
    // Get coords of Img(0,0,0)
    P.x=-(geo.sVoxelX/2-geo.dVoxelX/2)+geo.offOrigX[i];
    P.y=-(geo.sVoxelY/2-geo.dVoxelY/2)+geo.offOrigY[i];
    P.z=-(geo.sVoxelZ/2-geo.dVoxelZ/2)+geo.offOrigZ[i];
    
    // Get coors from next voxel in each direction
    Px.x=P.x+geo.dVoxelX;       Py.x=P.x;                Pz.x=P.x;
    Px.y=P.y;                   Py.y=P.y+geo.dVoxelY;    Pz.y=P.y;
    Px.z=P.z;                   Py.z=P.z;                Pz.z=P.z+geo.dVoxelZ;
    
    
    
   // Rotate image around X axis (this is equivalent of rotating the source and detector) RZ RY RZ
    eulerZYZT(geo,&P);
    eulerZYZT(geo,&Px);
    eulerZYZT(geo,&Py);
    eulerZYZT(geo,&Pz);
    
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
    
    
    Point3D source;
    source.x=0;
    source.y=-geo.offDetecU[i];
    source.z=-geo.offDetecV[i];
    
    rollPitchYawT(geo,i,&source);
    source.x=source.x-(geo.DSD[i]-geo.DSO[i]);
            
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
void CreateTextureParallel(float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream, bool alloc)
{
        //cudaArray Descriptor
#if IS_FOR_MATLAB_TIGRE
        const cudaExtent extent =make_cudaExtent(geo.nDetecV, geo.nDetecU, nangles);
#else
        const cudaExtent extent =make_cudaExtent(geo.nDetecU, geo.nDetecV, nangles);
#endif
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        //cuda Array
        if (alloc){
        cudaMalloc3DArray(&d_cuArrTex[0], &channelDesc, extent);
        cudaCheckErrors("Texture memory allocation fail");
        }
        cudaMemcpy3DParms copyParams = {0};
        
        
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)projectiondata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[0];
        copyParams.extent   = extent;
        copyParams.kind     = cudaMemcpyHostToDevice;
        cudaMemcpy3DAsync(&copyParams,stream[0+1]);
        cudaCheckErrors("Texture memory data copy fail");
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
        cudaCheckErrors("Texture object creation fail");
    
}