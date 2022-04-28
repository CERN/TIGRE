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
#include "voxel_backprojection.hpp"
#include "TIGRE_common.hpp"
#include <math.h>
#include "GpuIds.hpp"

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
    
    void CreateTexture(const GpuIds& gpuids,float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream, int nStreamDevice,bool allocate);


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

__constant__ Point3D projParamsArrayDev[6*PROJ_PER_KERNEL];  // Dev means it is on device

// We also need a corresponding array on the host side to be filled before each kernel call, then copied to the device (array in constant memory above)
// Point3D projParamsArrayHost[6*PROJ_PER_KERNEL];   // Host means it is host memory

// Now we also need to store sinAlpha and cosAlpha for each projection (two floats per projection)
__constant__ float projSinCosArrayDev[5*PROJ_PER_KERNEL];



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END RB, 10/31/2016: Add constant memory arrays to store parameters for all projections to be analyzed during a single kernel call
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
//
//      Function:       kernelPixelBackprojectionFDK
//
//      Description:    Main FDK backprojection kernel
//______________________________________________________________________________

__global__ void kernelPixelBackprojectionFDK(const Geometry geo, float* image,const int currProjSetNumber, const int totalNoOfProjections, cudaTextureObject_t tex)
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
        
        Point3D deltaX = projParamsArrayDev[6*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
        Point3D deltaY = projParamsArrayDev[6*projNumber+1];
        Point3D deltaZ = projParamsArrayDev[6*projNumber+2];
        Point3D xyzOrigin = projParamsArrayDev[6*projNumber+3];
        Point3D xyzOffset = projParamsArrayDev[6*projNumber+4];
        Point3D S = projParamsArrayDev[6*projNumber+5];
        
        float sinalpha = projSinCosArrayDev[5*projNumber];     // 2*projNumber because we have 2 float (sin or cos angle) values per projection
        float cosalpha = projSinCosArrayDev[5*projNumber+1];
        float COR = projSinCosArrayDev[5*projNumber+2];
        float DSD = projSinCosArrayDev[5*projNumber+3];
        float DSO = projSinCosArrayDev[5*projNumber+4];
        
        float auxCOR=COR/geo.dDetecU;
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
            P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y)-auxCOR;
            P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
            
            // This is the vector defining the line from the source to the Voxel
            float vectX,vectY,vectZ;
            vectX=(P.x -S.x);
            vectY=(P.y -S.y);
            vectZ=(P.z -S.z);
            
            // Get the coordinates in the detector UV where the mid point of the voxel is projected.
            float t=__fdividef(DSO-DSD-S.x,vectX);
            float y,z;
            y=vectY*t+S.y;
            z=vectZ*t+S.z;
            float u,v;
            u=y+(float)geo.nDetecU*0.5f;
            v=z+(float)geo.nDetecV*0.5f;
            
            float weigth;
            float realx,realy;
            realx=-(geo.sVoxelX-geo.dVoxelX)*0.5f  +indX*geo.dVoxelX   +xyzOffset.x;
            realy=-(geo.sVoxelY-geo.dVoxelY)*0.5f  +indY*geo.dVoxelY   +xyzOffset.y+COR;
            
            weigth=__fdividef(DSO+realy*sinalpha-realx*cosalpha,DSO);
            
            weigth=__frcp_rd(weigth*weigth);
            
            // Get Value in the computed (U,V) and multiply by the corresponding weigth.
            // indAlpha is the ABSOLUTE number of projection in the projection array (NOT the current number of projection set!)
            
#if IS_FOR_MATLAB_TIGRE
            voxelColumn[colIdx]+=tex3D<float>(tex, v, u ,indAlpha+0.5f)*weigth;
#else
            voxelColumn[colIdx]+=tex3D<float>(tex, u, v ,indAlpha+0.5f)*weigth;
#endif
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

int voxel_backprojection(float  *  projections, Geometry geo, float* result,float const * const alphas, int nalpha, const GpuIds& gpuids)
{
    // printf("voxel_backprojection(geo.nDetector = %d, %d)\n", geo.nDetecU, geo.nDetecV);
    // printf("geo.nVoxel    = %d, %d, %d\n", geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ);
    
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    cudaCheckErrors("Device query fail");
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Atb:Voxel_backprojection:GPUselect","There are no available device(s) that support CUDA\n");
    }

    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning thrown)
    // Check the available devices, and if they are the same
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Atb:Voxel_backprojection:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }

    int dev;
    // Split the CT problem
    unsigned int split_image;
    unsigned int split_projections;
    splitCTbackprojection(gpuids,geo,nalpha,&split_image,&split_projections);
    
    
    cudaCheckErrors("Error");
    //Pagelock memory for synchronous copy.
    // Lets try to make the host memory pinned:
    // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
    // empirical testing shows that when the image split is smaller than 1 (also implies the image is not very big), the time to
    // pin the memory is greater than the lost time in Synchronously launching the memcpys. This is only worth it when the image is too big.
    if (isHostRegisterSupported & (split_image>1 |deviceCount>1)){
        cudaHostRegister(result, (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*(size_t)sizeof(float),cudaHostRegisterPortable);
    }
//     if (isHostRegisterSupported ){
//         cudaHostRegister(projections, (size_t)geo.nDetecU*(size_t)geo.nDetecV*(size_t)nalpha*(size_t)sizeof(float),cudaHostRegisterPortable);
//     }
    cudaCheckErrors("Error pinning memory");
    
    
    // Create the arrays for the geometry. The main difference is that geo.offZ has been tuned for the
    // image slices. The rest of the Geometry is the same
    Geometry* geoArray=(Geometry*)malloc(split_image*deviceCount*sizeof(Geometry));
    createGeoArray(split_image*deviceCount,geo,geoArray,nalpha);
    
    // Now lest allocate all the image memory on the GPU, so we can use it later. If we have made our numbers correctly
    // in the previous section this should leave enough space for the textures.
    size_t num_bytes_img = (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geoArray[0].nVoxelZ* sizeof(float);
    float** dimage=(float**)malloc(deviceCount*sizeof(float*));
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMalloc((void**)&dimage[dev], num_bytes_img);
        cudaCheckErrors("cudaMalloc fail");
    }
    
    //If it is the first time, lets make sure our image is zeroed.
    int nStreamDevice=2;
    int nStreams=deviceCount*nStreamDevice;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));;
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        for (int i = 0; i < nStreamDevice; ++i){
            cudaStreamCreate(&stream[i+dev*nStreamDevice]);
            
        }
    }
    

     
    
    // Kernel auxiliary variables
    Point3D* projParamsArrayHost;
    cudaMallocHost((void**)&projParamsArrayHost,6*PROJ_PER_KERNEL*sizeof(Point3D));
    float* projSinCosArrayHost;
    cudaMallocHost((void**)&projSinCosArrayHost,5*PROJ_PER_KERNEL*sizeof(float));
    
    
    // Texture object variables
    cudaTextureObject_t *texProj;
    cudaArray **d_cuArrTex;
    texProj =(cudaTextureObject_t*)malloc(deviceCount*2*sizeof(cudaTextureObject_t));
    d_cuArrTex =(cudaArray**)malloc(deviceCount*2*sizeof(cudaArray*));
    
    // Auxiliary Host page-locked memory for fast and asycnornous memcpy.

    // Start with the main loop. The Projection data needs to be allocated and dealocated in the main loop
    // as due to the nature of cudaArrays, we can not reuse them. This should not be a problem for the fast execution
    // of the code, as repeated allocation and deallocation only happens when the projection data is very very big,
    // and therefore allcoation time should be negligible, fluctuation of other computations should mask the time.
    unsigned long long proj_linear_idx_start;
    unsigned int proj_split_overlap_number;
    unsigned int current_proj_split_size,current_proj_overlap_split_size;
    size_t num_bytes_img_curr;
    size_t img_linear_idx_start;
    float** partial_projection;
    size_t* proj_split_size;
    
    
    
    for(unsigned int img_slice=0;img_slice<split_image;img_slice++){
        // Initialize the memory if its the first time.
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaMemset(dimage[dev],0,num_bytes_img);
            cudaCheckErrors("memset fail");
        }
        
        for( unsigned int proj=0;proj<split_projections;proj++){
            
            
            // What is the size of the current chunk of proejctions we need in?
            current_proj_split_size=(nalpha+split_projections-1)/split_projections;
            // if its the last one its probably less
            current_proj_split_size=((proj+1)*current_proj_split_size<nalpha)?  current_proj_split_size:  nalpha-current_proj_split_size*proj;
            
            // We are going to split it in the same amount of kernels we need to execute.
            proj_split_overlap_number=(current_proj_split_size+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;
            
            // Create pointer to pointers of projections and precompute their location and size.
            if(!proj && !img_slice){
                partial_projection=(float**)malloc(proj_split_overlap_number*sizeof(float*));
                proj_split_size=(size_t*)malloc(proj_split_overlap_number*sizeof(size_t*));
            }
            for(unsigned int proj_block_split=0; proj_block_split<proj_split_overlap_number;proj_block_split++){
                // Crop the last one, as its likely its not completely divisible.
                // now lets split this for simultanoeus memcopy and compute.
                // We want to make sure that if we can, we run PROJ_PER_KERNEL projections, to maximize kernel acceleration
                // current_proj_overlap_split_size units = angles
                current_proj_overlap_split_size=max((current_proj_split_size+proj_split_overlap_number-1)/proj_split_overlap_number,PROJ_PER_KERNEL);
                current_proj_overlap_split_size=(proj_block_split<proj_split_overlap_number-1)?current_proj_overlap_split_size:current_proj_split_size-(proj_split_overlap_number-1)*current_proj_overlap_split_size;
                //Get the linear index where the current memory chunk starts.
                
                proj_linear_idx_start=(unsigned long long)((nalpha+split_projections-1)/split_projections)*(unsigned long long)proj*(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV;
                proj_linear_idx_start+=proj_block_split*max((current_proj_split_size+proj_split_overlap_number-1)/proj_split_overlap_number,PROJ_PER_KERNEL)*(unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV;
                //Store result
                proj_split_size[proj_block_split]=current_proj_overlap_split_size;
                partial_projection[proj_block_split]=&projections[proj_linear_idx_start];
                
            }                
            for(unsigned int proj_block_split=0; proj_block_split<proj_split_overlap_number;proj_block_split++){

                // Now get the projections on memory

                CreateTexture(gpuids,
                        partial_projection[proj_block_split],geo,
                        &d_cuArrTex[(proj_block_split%2)*deviceCount],
                        proj_split_size[proj_block_split],
                        &texProj   [(proj_block_split%2)*deviceCount],
                        stream, nStreamDevice,
                        (proj_block_split<2)&!proj&!img_slice);// Only allocate if its the first 2 calls
                
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(gpuids[dev]);
                    cudaStreamSynchronize(stream[dev*nStreamDevice+1]);
                 }
                               
                // Pin the next chunk of projection data, unpin the current one.
                for (dev = 0; dev < deviceCount; dev++){
                    //Safety:
                    // Depends on the amount of GPUs, the case where a image slice is zero hight can happen.
                    // Just break the loop if we reached that point
                    if(geoArray[img_slice*deviceCount+dev].nVoxelZ==0)
                        break;
                    
                    cudaSetDevice(gpuids[dev]);
                    
                    
                    
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
                    unsigned int noOfKernelCalls = (proj_split_size[proj_block_split]+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_KERNEL
                    for (unsigned int i=0; i<noOfKernelCalls; i++){
                        
                        // Now we need to generate and copy all data for PROJ_PER_KERNEL projections to constant memory so that our kernel can use it
                        unsigned int j;
                        for(j=0; j<PROJ_PER_KERNEL; j++){
                            
                            unsigned int currProjNumber_slice=i*PROJ_PER_KERNEL+j;
                            unsigned int currProjNumber_global=i*PROJ_PER_KERNEL+j                                                                          // index within kernel
                                                               +proj*(nalpha+split_projections-1)/split_projections                                          // index of the global projection split
                                                               +proj_block_split*max(current_proj_split_size/proj_split_overlap_number,PROJ_PER_KERNEL); // indexof overlap current split
                            
                            if(currProjNumber_slice>=proj_split_size[proj_block_split])
                                break;  // Exit the loop. Even when we leave the param arrays only partially filled, this is OK, since the kernel will check bounds anyway.
                            if(currProjNumber_global>=nalpha)
                                break;  // Exit the loop. Even when we leave the param arrays only partially filled, this is OK, since the kernel will check bounds anyway.
                            
                            Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, /*offDetec,*/source;
                            float sinalpha,cosalpha;
                            
                            geoArray[img_slice*deviceCount+dev].alpha=-alphas[currProjNumber_global*3];//we got 3 angles now.
                            geoArray[img_slice*deviceCount+dev].theta=-alphas[currProjNumber_global*3+1];
                            geoArray[img_slice*deviceCount+dev].psi  =-alphas[currProjNumber_global*3+2];
                            
//                             mexPrintf("%u %f \n",i,geoArray[img_slice*deviceCount+dev].alpha);
//                             mexPrintf("%u \n",currProjNumber_global);
                            
                            sinalpha=sin(geoArray[img_slice*deviceCount+dev].alpha);
                            cosalpha=cos(geoArray[img_slice*deviceCount+dev].alpha);
                            
                            projSinCosArrayHost[5*j]=sinalpha;  // 2*j because we have 2 float (sin or cos angle) values per projection
                            projSinCosArrayHost[5*j+1]=cosalpha;
                            projSinCosArrayHost[5*j+2]=geo.COR[currProjNumber_global];
                            projSinCosArrayHost[5*j+3]=geo.DSD[currProjNumber_global];
                            projSinCosArrayHost[5*j+4]=geo.DSO[currProjNumber_global];
                            
                            computeDeltasCube(geoArray[img_slice*deviceCount+dev],currProjNumber_global,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
                            
                            offOrig.x=geo.offOrigX[currProjNumber_global];
                            offOrig.y=geo.offOrigY[currProjNumber_global];
                            offOrig.z=geoArray[img_slice*deviceCount+dev].offOrigZ[currProjNumber_global];
                            
                            projParamsArrayHost[6*j]=deltaX;		// 6*j because we have 6 Point3D values per projection
                            projParamsArrayHost[6*j+1]=deltaY;
                            projParamsArrayHost[6*j+2]=deltaZ;
                            projParamsArrayHost[6*j+3]=xyzOrigin;
                            projParamsArrayHost[6*j+4]=offOrig;
                            projParamsArrayHost[6*j+5]=source;
                        }   // END for (preparing params for kernel call)
                        
                        // Copy the prepared parameter arrays to constant memory to make it available for the kernel
                        cudaMemcpyToSymbolAsync(projSinCosArrayDev, projSinCosArrayHost, sizeof(float)*5*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[dev*nStreamDevice]);
                        cudaMemcpyToSymbolAsync(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*6*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[dev*nStreamDevice]);
                        cudaStreamSynchronize(stream[dev*nStreamDevice]);
                        
                        kernelPixelBackprojectionFDK<<<grid,block,0,stream[dev*nStreamDevice]>>>(geoArray[img_slice*deviceCount+dev],dimage[dev],i,proj_split_size[proj_block_split],texProj[(proj_block_split%2)*deviceCount+dev]);
                    }  // END for
                    //////////////////////////////////////////////////////////////////////////////////////
                    // END RB code, Main reconstruction loop: go through projections (rotation angles) and backproject
                    //////////////////////////////////////////////////////////////////////////////////////
                }// END for deviceCount
            } // END sub-split of current projection chunk
            
            for (dev = 0; dev < deviceCount; dev++){
                cudaSetDevice(gpuids[dev]);
                cudaDeviceSynchronize();
            }
            
        } // END projection splits
        
       
        // Now we need to take the image out of the GPU
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            // We do not need to sycnronize because the array dealocators already do.
            num_bytes_img_curr=(size_t)geoArray[img_slice*deviceCount+dev].nVoxelX*(size_t)geoArray[img_slice*deviceCount+dev].nVoxelY*(size_t)geoArray[img_slice*deviceCount+dev].nVoxelZ*sizeof(float);
            img_linear_idx_start=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geoArray[0].nVoxelZ*(size_t)(img_slice*deviceCount+dev);
            cudaMemcpyAsync(&result[img_linear_idx_start], dimage[dev], num_bytes_img_curr, cudaMemcpyDeviceToHost,stream[dev*nStreamDevice+1]);
        }
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
            cudaCheckErrors("Main loop fail");
        }
        
    } // end image splits

    ///////// Cleaning:
    
    
    bool two_buffers_used=((((nalpha+split_projections-1)/split_projections)+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL)>1;
    for(unsigned int i=0; i<2;i++){ // 2 buffers (if needed, maybe only 1)
        if (!two_buffers_used && i==1)
            break;
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDestroyTextureObject(texProj[i*deviceCount+dev]);
            cudaFreeArray(d_cuArrTex[i*deviceCount+dev]);
        }
    }
    cudaCheckErrors("cudadestroy textures result fail");
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaFree(dimage[dev]);
    }
    cudaFreeHost(projSinCosArrayHost);
    cudaFreeHost(projParamsArrayHost);
    free(partial_projection);
    free(proj_split_size);
    
    freeGeoArray(split_image*deviceCount,geoArray);
    
    if (isHostRegisterSupported & (split_image>1 |deviceCount>1)){
        cudaHostUnregister(result);
    }
//     if (isHostRegisterSupported){
//         cudaHostUnregister(projections);
//     }
//     
    
    for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]);
    
    cudaCheckErrors("cudaFree fail");
    
    cudaDeviceReset(); // For the Nvidia Visual Profiler
    return 0;
    
}  // END voxel_backprojection
//

void splitCTbackprojection(const GpuIds& gpuids, Geometry geo,int nalpha, unsigned int* split_image, unsigned int * split_projections){
    
    
    // We don't know if the devices are being used. lets check that. and only use the amount of memory we need.
    
    size_t mem_GPU_global;
    checkFreeMemory(gpuids, &mem_GPU_global);

    const int deviceCount = gpuids.GetLength();
        
    // Compute how much memory each of the relevant memory pieces need
    size_t mem_image=       (unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
    size_t mem_proj=        (unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV*sizeof(float);
    
    
    
    
    // Does everything fit in the GPU?
    
    if(mem_image/deviceCount+mem_proj*PROJ_PER_KERNEL*2<mem_GPU_global){
        // We only need to split if we have extra GPUs
        *split_image=1;
        *split_projections=1;
    }
    // We know we need to split, but:
    // Does all the image fit in the GPU, with some slack for a stack of projections??
    else
    {
        // As we can overlap memcpys from H2D of the projections, we should then minimize the amount of image splits.
        // Lets assume to start with that we only need 1 stack of PROJ_PER_KERNEL projections. The rest is for the image.
        size_t mem_free=mem_GPU_global-2*mem_proj*PROJ_PER_KERNEL;
        
        *split_image=(mem_image/deviceCount+mem_free-1)/mem_free;
        // Now knowing how many splits we have for images, we can recompute how many slices of projections actually
        // fit on the GPU. Must be more than 0 obviously.
        
        mem_free=mem_GPU_global-(mem_image/deviceCount)/(*split_image); // NOTE: There is some rounding error, but its in the order of bytes, and we have 5% of GPU free jsut in case. We are safe
        
        
        *split_projections=(mem_proj*PROJ_PER_KERNEL*2+mem_free-1)/mem_free;
        
    }
}


void CreateTexture(const GpuIds& gpuids, float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream,int nStreamDevice,bool allocate){
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
#if IS_FOR_MATLAB_TIGRE
    const cudaExtent extent =make_cudaExtent(geo.nDetecV, geo.nDetecU, nangles);
#else
    const cudaExtent extent =make_cudaExtent(geo.nDetecU, geo.nDetecV, nangles);
#endif
    const unsigned int num_devices = gpuids.GetLength();
    if (allocate){
        for (unsigned int dev = 0; dev < num_devices; dev++){
            cudaSetDevice(gpuids[dev]);
            
            //cudaArray Descriptor
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            //cuda Array
            cudaMalloc3DArray(&d_cuArrTex[dev], &channelDesc, extent);
            
        }
    }
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemcpy3DParms copyParams = {0};
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)projectiondata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[dev];
        copyParams.extent   = extent;
        copyParams.kind     = cudaMemcpyHostToDevice;
        cudaMemcpy3DAsync(&copyParams,stream[dev*nStreamDevice+1]);
    }

    //Array creation End
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArrTex[dev];
        cudaTextureDesc     texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = false;
        texDescr.filterMode = cudaFilterModeLinear;
        texDescr.addressMode[0] = cudaAddressModeBorder;
        texDescr.addressMode[1] = cudaAddressModeBorder;
        texDescr.addressMode[2] = cudaAddressModeBorder;
        texDescr.readMode = cudaReadModeElementType;
        cudaCreateTextureObject(&texImage[dev], &texRes, &texDescr, NULL);
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
//______________________________________________________________________________
//
//      Function:       computeDeltasCube
//
//      Description:    Computes relative increments for each projection (volume rotation).
//						Increments get passed to the backprojection kernel.
//______________________________________________________________________________

void computeDeltasCube(Geometry geo,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S)
{
    
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
    //Done for P, now source
    Point3D source;
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

void eulerZYZT(Geometry geo, Point3D* point){
    
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x = auxPoint.x*(cos(geo.psi)*cos(geo.theta)*cos(geo.alpha)-sin(geo.psi)*sin(geo.alpha))
    +auxPoint.y*(-cos(geo.psi)*cos(geo.theta)*sin(geo.alpha)-sin(geo.psi)*cos(geo.alpha))
    +auxPoint.z*cos(geo.psi)*sin(geo.theta);
    point->y = auxPoint.x*(sin(geo.psi)*cos(geo.theta)*cos(geo.alpha)+cos(geo.psi)*sin(geo.alpha))
    +auxPoint.y*(-sin(geo.psi)*cos(geo.theta)*sin(geo.alpha)+cos(geo.psi)*cos(geo.alpha))
    +auxPoint.z*sin(geo.psi)*sin(geo.theta);
    point->z =-auxPoint.x*sin(geo.theta)*cos(geo.alpha)
    +auxPoint.y*sin(geo.theta)*sin(geo.alpha)
    +auxPoint.z*cos(geo.theta);
}
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
void checkFreeMemory(const GpuIds& gpuids,size_t *mem_GPU_global){
    size_t memfree;
    size_t memtotal;
    const int deviceCount = gpuids.GetLength();
    
    for (int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) *mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            mexErrMsgIdAndTxt("tvDenoise:tvdenoising:GPU","One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n Free the GPU to run TIGRE\n");
        }
        cudaCheckErrors("Check mem error");
        
        *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
    }
    *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
    
    //*mem_GPU_global= insert your known number here, in bytes.
}

