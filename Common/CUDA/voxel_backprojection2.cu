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
#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include "page_locking.hpp"
#include "launch_shapes.hpp"
#include <math.h>
#include "GpuIds.hpp"

// https://stackoverflow.com/questions/16282136/is-there-a-cuda-equivalent-of-perror
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                tigreCudaError("CBCT:CUDA:Atb", msg, cudaGetErrorString(__err));\
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
void CreateTexture2(const GpuIds& gpuids, float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream,int nStreamDevice,bool allocate, cudaMemcpyKind kind);
static void matchedKernelParams(const Geometry& geo, int i, Point3D* params);  // see its definition

__global__ void matrixConstantMultiply(const Geometry geo,float* image,float constant){
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    for(; idx<(size_t)geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ; idx+=gridDim.x*blockDim.x) {
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

// Now we also need to store sinAlpha and cosAlpha for each projection (two floats per projection)
__constant__ float projSinCosArray2Dev[5*PROJ_PER_KERNEL];

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

// VPT: voxels per thread along z, VOXELS_PER_THREAD or fewer on thin images (tigre::voxelColumnLength).
template<int VPT>
__global__ void kernelPixelBackprojection(const Geometry geo, float* image,const int currProjSetNumber, const int totalNoOfProjections, cudaTextureObject_t tex)
{
    
    unsigned long long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long long indX = blockIdx.x * blockDim.x + threadIdx.x;
    // unsigned long startIndZ = blockIdx.z * blockDim.z + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    unsigned long long startIndZ = blockIdx.z * VPT + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    //Make sure we don't go out of bounds
    if (indX>=geo.nVoxelX || indY>=geo.nVoxelY || startIndZ>=geo.nVoxelZ)
        return;
    
    // We'll keep a local auxiliary array of values of a column of voxels that this thread will update
    // Zero-initialised: voxels of the column beyond the image are computed too (no per-voxel check, see below), never stored.
    float voxelColumn[VPT] = {0};
    
    // First we need to copy the curent 3D volume values from the column to our auxiliary array so that we can then
    // work on them (update them by computing values from multiple projections) locally - avoiding main memory reads/writes
    
    unsigned long colIdx;
    // The voxel's indices as floats, converted once per thread (they were converted for every voxel and projection).
    const float fIndX=(float)indX, fIndY=(float)indY, fStartIndZ=(float)startIndZ;
#pragma unroll
    for(colIdx=0; colIdx<VPT; colIdx++)
    {
        unsigned long long indZ = startIndZ + colIdx;
        // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
        // be trying to copy the out of bounds values back to the 3D volume anyway (bounds checks will be done in the final loop where the updated values go back to the main volume)
        if(indZ>=geo.nVoxelZ)
            break;   // break the loop.
        
        unsigned long long idx =indZ*(unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY+indY*(unsigned long long)geo.nVoxelX + indX;
        voxelColumn[colIdx] = image[idx];   // Read the current volume value that we'll update by computing values from MULTIPLE projections (not just one)
        // We'll be updating the local (register) variable, avoiding reads/writes from the slow main memory.
    }  // END copy 3D volume voxels to local array
    
    // Now iterate through projections
#pragma unroll
    for(unsigned long projNumber=0; projNumber<PROJ_PER_KERNEL; projNumber++)
    {
        // Get the current parameters from parameter arrays in constant memory.
        unsigned long indAlpha = currProjSetNumber*PROJ_PER_KERNEL+projNumber;  // This is the ABSOLUTE projection number in the projection array
        
        // Our currImageVal will be updated by hovewer many projections we had left in the "remainder" - that's OK.
        if(indAlpha>=totalNoOfProjections)
            break;
        
        // The projection's parameters, folded on the host (see the main loop): the source-to-voxel vector is
        // c0 + i*deltaX + j*deltaY + k*deltaZ for voxel (i,j,k), in the scaled coordinate system rotated with the
        // projection (detector rotation, offsets and COR included).
        Point3D deltaX = projParamsArray2Dev[7*projNumber];  // 7*projNumber because we have 7 Point3D values per projection
        Point3D deltaY = projParamsArray2Dev[7*projNumber+1];
        Point3D deltaZ = projParamsArray2Dev[7*projNumber+2];
        Point3D c0     = projParamsArray2Dev[7*projNumber+3];
        Point3D tuv    = projParamsArray2Dev[7*projNumber+4];  // (DSO-DSD-S.x, S.y+nDetecU/2, S.z+nDetecV/2)
        Point3D wpar   = projParamsArray2Dev[7*projNumber+5];  // (1/d_n, dDetecU^2, dDetecV^2) for the weight
        // "XYZ" in the scaled coordinate system, rotated with the projection angle, is affine in the voxel's indices,
        // and so is the vector from the source to the voxel: compute it for the first voxel of the column, then step
        // it by deltaZ (at most VPT-1 steps). It was rebuilt for every voxel from the indices.
        float vectX=c0.x+fIndX*deltaX.x+fIndY*deltaY.x+fStartIndZ*deltaZ.x;
        float vectY=c0.y+fIndX*deltaX.y+fIndY*deltaY.y+fStartIndZ*deltaZ.y;
        float vectZ=c0.z+fIndX*deltaX.z+fIndY*deltaY.z+fStartIndZ*deltaZ.z;
        const float tNum=tuv.x, u0=tuv.y, v0=tuv.z;
        // indAlpha is the ABSOLUTE number of projection in the projection array (NOT the current number of projection set!)
        const float fAlpha=(float)indAlpha+0.5f;
        
        // The weight (see below). In this frame the detector is the plane x=DSO-DSD, whatever the detector's rotation, and
        // only y and z are scaled (to detector pixels): -tNum is the distance from the source to the detector plane, and
        // the source-to-voxel distance is |(vectX, vectY*dDetecU, vectZ*dDetecV)|.
        const float invDn=wpar.x, dU2=wpar.y, dV2=wpar.z;
        
        // Now iterate through Z in our voxel column FOR A GIVEN PROJECTION
#pragma unroll
        // No bounds check per voxel: the voxels of a column that lie beyond the image (in the last block in z) are
        // computed too, harmlessly (a texture read), and never stored (the final loop checks).
        for(colIdx=0; colIdx<VPT; colIdx++)
        {
            // Get the coordinates in the detector UV where the mid point of the voxel is projected.
            float t=__fdividef(tNum,vectX);
            float u=vectY*t+u0;
            float v=vectZ*t+v0;
#if IS_FOR_MATLAB_TIGRE
            float sample=tex3D<float>(tex, v, u ,fAlpha);
#else
            float sample=tex3D<float>(tex, u, v ,fAlpha);
#endif
            // Weight L^3/(d_n*|SV|^2): L the distance from the source to where the ray through the voxel meets the
            // detector, |SV| from the source to the voxel, d_n from the source to the detector plane (DSD if the
            // detector is not rotated). The ray meets the detector at t times the source-to-voxel vector, so
            // L=t*|SV| and the weight is t^3*|SV|/d_n. Before, L was rebuilt in world coordinates from (u,v) as if the
            // detector were not rotated (and the trajectory a circle about z), half a pixel off: off by ~2e-4 for an
            // unrotated detector and by several % for a rotated one, and about 45 instructions per voxel.
            float lsq=vectX*vectX+dU2*vectY*vectY+dV2*vectZ*vectZ;
            float weight=t*t*t*(lsq*rsqrtf(lsq))*invDn;
            voxelColumn[colIdx]+=sample* weight;
            vectX+=deltaZ.x; vectY+=deltaZ.y; vectZ+=deltaZ.z;
        }  // END iterating through column of voxels
        
    }  // END iterating through multiple projections
    
    // And finally copy the updated local voxelColumn array back to our 3D volume (main memory)
#pragma unroll
    for(colIdx=0; colIdx<VPT; colIdx++)
    {
        unsigned long long indZ = startIndZ + colIdx;
        // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
        // be trying to copy the out of bounds values back to the 3D volume anyway (bounds checks will be done in the final loop where the values go to the main volume)
        if(indZ>=geo.nVoxelZ)
            break;   // break the loop.
        
        unsigned long long idx =indZ*(unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY+indY*(unsigned long long)geo.nVoxelX + indX;
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

int voxel_backprojection2(float * projections, Geometry geo, float* result,float const * const alphas, int nalpha, const GpuIds& gpuids, bool gpuData){
    
    
    
    
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    cudaCheckErrors("Device query fail");
    // If a CUDA error raises (Python build) before the end of this function,
    // onError releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    if (gpuData && deviceCount != 1) {
        mexErrMsgIdAndTxt("Atb:GPUdata", "GPU arrays can only be backprojected on one GPU, the one they are on.");
    }
    std::vector<int> devIds(deviceCount);  // copied: the release actions must not hold on to gpuids
    for (int d = 0; d < deviceCount; d++) { devIds[d] = gpuids[d]; onError.addDevice(devIds[d]); }
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Atb:Voxel_backprojection:GPUselect","There are no available device(s) that support CUDA\n");
    }
    
    
    // CODE assumes
    // 1.-All available devices are usable by this code
    // 2.-All available devices are equal, they are the same machine (warning thrown)
    // Check the available devices, and if they are the same
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Atb:Voxel_backprojection2:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }

    int dev;

    
    // Split the CT problem
    unsigned int split_image;
    unsigned int split_projections;
    bool one_upload;
    splitCTbackprojection(gpuids,geo,nalpha,&split_image,&split_projections,&one_upload);
    
    
    // Create the arrays for the geometry. The main difference is that geo.offZ has been tuned for the
    // image slices. The rest of the Geometry is the same
    Geometry* geoArray=(Geometry*)malloc(split_image*deviceCount*sizeof(Geometry));
    createGeoArray(split_image*deviceCount,geo,geoArray,nalpha);
    onError.add([geoArray, split_image, deviceCount]{ freeGeoArray(split_image*deviceCount, geoArray); });
    
    // Now lest allocate all the image memory on the GPU, so we can use it later. If we have made our numbers correctly
    // in the previous section this should leave enough space for the textures.
    size_t num_bytes_img = (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geoArray[0].nVoxelZ* sizeof(float);
    float** dimage=(float**)malloc(deviceCount*sizeof(float*));
    memset(dimage, 0, deviceCount*sizeof(float*));
    // GPU data, whole image in one piece: backproject straight into the caller's
    // array instead of into a buffer that is then copied there.
    const bool directOut = gpuData && split_image == 1 && deviceCount == 1;
    onError.add([dimage, deviceCount, devIds, result]{
        for (int d = 0; d < deviceCount; d++) {
            cudaSetDevice(devIds[d]);
            if (dimage[d] != result) cudaFree(dimage[d]);
        }
        free(dimage);
    });
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        if (directOut)
            dimage[dev] = result;
        else
            cudaMalloc((void**)&dimage[dev], num_bytes_img);
        cudaCheckErrors("cudaMalloc fail");
    }
        
    
    //Pagelock memory for synchronous copy.
    // Lets try to make the host memory pinned:
    // We laredy queried the GPU and assuemd they are the same, thus should have the same attributes.
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,gpuids[0]);
#endif
    if (gpuData) isHostRegisterSupported = 0;  // GPU memory: nothing to pin
    // empirical testing shows that when the image split is smaller than 1 (also implies the image is not very big), the time to
    // pin the memory is greater than the lost time in Synchronously launching the memcpys. This is only worth it when the image is too big.
    // tigre::hostRegister: an array that is already page-locked (e.g. pinned by PyTorch) is used as it is and not unregistered here.
    bool result_registered=false, projections_registered=false;
    if (tigre::worthPageLocking(isHostRegisterSupported, gpuData, split_image, deviceCount)){  // as FDK: also with several GPUs
        result_registered=tigre::hostRegister(result, (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*(size_t)sizeof(float));
        if (result_registered) onError.add([result]{ cudaHostUnregister(result); });
    }
    // The projections cross more than once only with several GPUs (each gets all of them) or a split image (they are
    // sent again for each piece). Sent once, their chunks are copied while the GPU works on the previous one anyway.
    if (tigre::worthPageLocking(isHostRegisterSupported, gpuData, split_image, deviceCount)){
        projections_registered=tigre::hostRegister(projections, (size_t)geo.nDetecU*(size_t)geo.nDetecV*(size_t)nalpha*(size_t)sizeof(float));
        // The caller's array: left registered, the next call with it would fail.
        if (projections_registered) onError.add([projections]{ cudaHostUnregister(projections); });
    }
    cudaCheckErrors("Error pinning memory");

    
    
    

    //If it is the first time, lets make sure our image is zeroed.
    int nStreamDevice=2;
    int nStreams=deviceCount*nStreamDevice;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    memset(stream, 0, nStreams*sizeof(cudaStream_t));
    onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        for (int i = 0; i < nStreamDevice; ++i){
            cudaStreamCreate(&stream[i+dev*nStreamDevice]);
            
        }
    }
    
    // Kernel auxiliary variables
    Point3D* projParamsArray2Host = 0;
    cudaMallocHost((void**)&projParamsArray2Host,7*PROJ_PER_KERNEL*sizeof(Point3D));
    onError.add([projParamsArray2Host]{ cudaFreeHost(projParamsArray2Host); });
    float* projSinCosArray2Host = 0;
    cudaMallocHost((void**)&projSinCosArray2Host,5*PROJ_PER_KERNEL*sizeof(float));
    onError.add([projSinCosArray2Host]{ cudaFreeHost(projSinCosArray2Host); });
    
    // Texture object variables
    cudaTextureObject_t *texProj;
    cudaArray **d_cuArrTex;
    texProj =(cudaTextureObject_t*)malloc(deviceCount*2*sizeof(cudaTextureObject_t));
    // Zeroed so CreateTexture can tell an empty slot from one holding a texture
    // object it must destroy first.
    memset(texProj, 0, deviceCount*2*sizeof(cudaTextureObject_t));
    d_cuArrTex =(cudaArray**)malloc(deviceCount*2*sizeof(cudaArray*));
    memset(d_cuArrTex, 0, deviceCount*2*sizeof(cudaArray*));
    onError.add([texProj, d_cuArrTex, deviceCount, devIds]{
        for (int s = 0; s < 2*deviceCount; s++) {
            cudaSetDevice(devIds[s % deviceCount]);
            if (texProj[s]) cudaDestroyTextureObject(texProj[s]);
            if (d_cuArrTex[s]) cudaFreeArray(d_cuArrTex[s]);
        }
        free(texProj); free(d_cuArrTex);
    });
    
    
    
    unsigned int proj_split_overlap_number;
    // Start with the main loop. The Projection data needs to be allocated and dealocated in the main loop
    // as due to the nature of cudaArrays, we can not reuse them. This should not be a problem for the fast execution
    // of the code, as repeated allocation and deallocation only happens when the projection data is very very big,
    // and therefore allcoation time should be negligible, fluctuation of other computations should mask the time.
    unsigned long long proj_linear_idx_start;
    unsigned int current_proj_split_size,current_proj_overlap_split_size;
    size_t num_bytes_img_curr;
    size_t img_linear_idx_start;
    float** partial_projection = 0;
    size_t* proj_split_size = 0;
    
    for(unsigned int img_slice=0;img_slice<split_image;img_slice++){
//
        // Initialize the memory if its the first time.
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaMemsetAsync(dimage[dev],0,num_bytes_img,stream[dev*nStreamDevice]); // on the kernels' stream: ordered before them
            cudaCheckErrors("memset fail");
        }
        
        for( unsigned int proj=0;proj<split_projections;proj++){
            
            
            // What is the size of the current chunk of proejctions we need in?
            current_proj_split_size=(nalpha+split_projections-1)/split_projections;
            // if its the last one its probably less
            current_proj_split_size=((proj+1)*current_proj_split_size<nalpha)?  current_proj_split_size:  nalpha-current_proj_split_size*proj;
            
            // We are going to split it in the same amount of kernels we need to execute.
            // (all of them in one go on thin detectors: see splitCTbackprojection)
            proj_split_overlap_number=one_upload ? 1 : (current_proj_split_size+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;
            
            // Create pointer to pointers of projections and precompute their location and size.
            if(!proj && !img_slice){
                partial_projection=(float**)malloc(current_proj_split_size*sizeof(float*));
                proj_split_size=(size_t*)malloc(current_proj_split_size*sizeof(size_t*));
                onError.add([partial_projection, proj_split_size]{ free(partial_projection); free(proj_split_size); });
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

                CreateTexture2(gpuids,
                        partial_projection[proj_block_split],geo,
                        &d_cuArrTex[(proj_block_split%2)*deviceCount],
                        proj_split_size[proj_block_split],
                        &texProj   [(proj_block_split%2)*deviceCount],
                        stream, nStreamDevice,
                        (proj_block_split<2)&!proj&!img_slice,// Only allocate if its the first 2 calls
                        (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice));
                // Stop here if the texture could not be made (e.g. out of memory),
                // rather than run the kernels on a broken texture.
                cudaCheckErrors("Texture creation fail");
                
                for (dev = 0; dev < deviceCount; dev++){
                    cudaSetDevice(gpuids[dev]);
                    cudaStreamSynchronize(stream[dev*nStreamDevice+1]);
                 }

                for (dev = 0; dev < deviceCount; dev++){
                    //Safety:
                    // Depends on the amount of GPUs, the case where a image slice is zero hight can happen.
                    // Just break the loop if we reached that point
                    if(geoArray[img_slice*deviceCount+dev].nVoxelZ==0)
                        break;
                    
                    cudaSetDevice(gpuids[dev]);
                    
                    
                    
                    int divx,divy,divz;
                    // RB: Use the optimal (in their tests) block size from paper by Zinsser and Keck (16 in x and 32 in y).
                    // I tried different sizes and shapes of blocks (tiles), but it does not appear to significantly affect throughput, so
                    // let's stick with the values from Zinsser and Keck.
                    divx=16;
                    divy=32;
                    divz=tigre::voxelColumnLength(geoArray[img_slice*deviceCount+dev].nVoxelZ, VOXELS_PER_THREAD);  // We now only have 32 x 16 threads per block (flat tile, see below), BUT each thread works on a Z column of divz voxels (VOXELS_PER_THREAD, or fewer on thin images such as fan beam), so we effectively need fewer blocks!
                    
                    
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
                            
                            float sinalpha,cosalpha;
                            
                            geoArray[img_slice*deviceCount+dev].alpha=-alphas[currProjNumber_global*3];//we got 3 angles now.
                            geoArray[img_slice*deviceCount+dev].theta=-alphas[currProjNumber_global*3+1];
                            geoArray[img_slice*deviceCount+dev].psi  =-alphas[currProjNumber_global*3+2];
                            
                            sinalpha=sin(geoArray[img_slice*deviceCount+dev].alpha);
                            cosalpha=cos(geoArray[img_slice*deviceCount+dev].alpha);
                            
                            projSinCosArray2Host[5*j]=sinalpha;  // 2*j because we have 2 float (sin or cos angle) values per projection
                            projSinCosArray2Host[5*j+1]=cosalpha;
                            projSinCosArray2Host[5*j+2]=geo.COR[currProjNumber_global];
                            projSinCosArray2Host[5*j+3]=geo.DSD[currProjNumber_global];
                            projSinCosArray2Host[5*j+4]=geo.DSO[currProjNumber_global];
                            
                            matchedKernelParams(geoArray[img_slice*deviceCount+dev],currProjNumber_global,&projParamsArray2Host[7*j]);  // 7 Point3D per projection
                            
                        }   // END for (preparing params for kernel call)
                        
                        // Copy the prepared parameter arrays to constant memory to make it available for the kernel
                        cudaMemcpyToSymbolAsync(projSinCosArray2Dev, projSinCosArray2Host, sizeof(float)*5*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[dev*nStreamDevice]);
                        cudaMemcpyToSymbolAsync(projParamsArray2Dev, projParamsArray2Host, sizeof(Point3D)*7*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[dev*nStreamDevice]);
                        cudaStreamSynchronize(stream[dev*nStreamDevice]);
                        switch (divz) {  // one compiled kernel per column length
                            case 1: kernelPixelBackprojection<1><<<grid,block,0,stream[dev*nStreamDevice]>>>(geoArray[img_slice*deviceCount+dev],dimage[dev],i,proj_split_size[proj_block_split],texProj[(proj_block_split%2)*deviceCount+dev]); break;
                            case 2: kernelPixelBackprojection<2><<<grid,block,0,stream[dev*nStreamDevice]>>>(geoArray[img_slice*deviceCount+dev],dimage[dev],i,proj_split_size[proj_block_split],texProj[(proj_block_split%2)*deviceCount+dev]); break;
                            case 4: kernelPixelBackprojection<4><<<grid,block,0,stream[dev*nStreamDevice]>>>(geoArray[img_slice*deviceCount+dev],dimage[dev],i,proj_split_size[proj_block_split],texProj[(proj_block_split%2)*deviceCount+dev]); break;
                            default: kernelPixelBackprojection<VOXELS_PER_THREAD><<<grid,block,0,stream[dev*nStreamDevice]>>>(geoArray[img_slice*deviceCount+dev],dimage[dev],i,proj_split_size[proj_block_split],texProj[(proj_block_split%2)*deviceCount+dev]); break;
                        }
                        
                    }  // END for
                    //////////////////////////////////////////////////////////////////////////////////////
                    // END RB code, Main reconstruction loop: go through projections (rotation angles) and backproject
                    //////////////////////////////////////////////////////////////////////////////////////
                }
            } // END sub-split of current projection chunk

        } // END projection splits
        
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            matrixConstantMultiply<<<60,MAXTREADS,0,stream[dev*nStreamDevice]>>>(  geoArray[img_slice*deviceCount+dev],dimage[dev],geo.dVoxelX*geo.dVoxelY*geo.dVoxelZ/(geo.dDetecU*geo.dDetecV));
        }

        // Now we need to take the image out of the GPU
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaStreamSynchronize(stream[dev*nStreamDevice]);
            
            num_bytes_img_curr=(size_t)geoArray[img_slice*deviceCount+dev].nVoxelX*(size_t)geoArray[img_slice*deviceCount+dev].nVoxelY*(size_t)geoArray[img_slice*deviceCount+dev].nVoxelZ*sizeof(float);
            img_linear_idx_start=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geoArray[0].nVoxelZ*(size_t)(img_slice*deviceCount+dev);
            if (!directOut)
                cudaMemcpyAsync(&result[img_linear_idx_start], dimage[dev], num_bytes_img_curr, (gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost),stream[dev*nStreamDevice+1]);
        }
        // The copy is asynchronous (to GPU memory, or to host memory pinned
        // above when the image is split or several GPUs are used): finish it before the next image
        // piece clears dimage. Without this, the cudaMemset of the next piece
        // raced the copy and parts of the result came out as zeros.
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDeviceSynchronize();
        }
    } // end image splits
    
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaDeviceSynchronize();
    }  
    cudaCheckErrors("Main loop fail");
    onError.dismiss();  // from here on the function releases everything itself
    
    // Clean the GPU
    bool two_buffers_used=!one_upload && ((((nalpha+split_projections-1)/split_projections)+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL)>1;
    for(unsigned int i=0; i<2;i++){ // 2 buffers (if needed, maybe only 1)
        if (!two_buffers_used && i==1)
            break;        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(gpuids[dev]);
            cudaDestroyTextureObject(texProj[i*deviceCount+dev]);
            cudaFreeArray(d_cuArrTex[i*deviceCount+dev]);
        }
    }
    free(d_cuArrTex);
    free(texProj);

    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        if (!directOut)
            cudaFree(dimage[dev]);
    }
    free(dimage);
    
    cudaFreeHost(projSinCosArray2Host);
    cudaFreeHost(projParamsArray2Host);
    free(partial_projection);
    free(proj_split_size);
    
    freeGeoArray(split_image*deviceCount,geoArray);
#ifndef NO_PINNED_MEMORY     
    if (result_registered){
        cudaHostUnregister(result);
    }
    if (projections_registered){
        cudaHostUnregister(projections);
    }
#endif 
    for (int i = 0; i < nStreams; ++i)
        cudaStreamDestroy(stream[i]);
    
    cudaCheckErrors("cudaFree fail");
    
//     cudaDeviceReset(); // For the Nvidia Visual Profiler
    return 0;
    
}  // END voxel_backprojection





void CreateTexture2(const GpuIds& gpuids, float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream,int nStreamDevice,bool allocate, cudaMemcpyKind kind){
    //size_t size_image=geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ;
    int num_devices = gpuids.GetLength();
#if IS_FOR_MATLAB_TIGRE
    const cudaExtent extent =make_cudaExtent(geo.nDetecV, geo.nDetecU, nangles);
#else
    const cudaExtent extent =make_cudaExtent(geo.nDetecU, geo.nDetecV, nangles);
#endif
    if (allocate){
        for (unsigned int dev = 0; dev < num_devices; dev++){
            cudaSetDevice(gpuids[dev]);
            
            //cudaArray Descriptor
            cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
            //cuda Array
            tigre::checkMalloc(cudaMalloc3DArray(&d_cuArrTex[dev], &channelDesc, extent), "Atb:voxel_backprojection2:GPUmemory", "Backprojection", (size_t)extent.width*extent.height*extent.depth*sizeof(float));
            
        }
    }
    for (unsigned int dev = 0; dev < num_devices; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemcpy3DParms copyParams = {0};
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)projectiondata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[dev];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
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
        // Each projection block (or image split) creates a new texture object in this
        // slot. Destroy the one it replaces, otherwise every call leaks one per block
        // until cudaCreateTextureObject fails with "out of memory" (#617). The kernels
        // that used it have finished: its array has just been overwritten above.
        if (texImage[dev]) cudaDestroyTextureObject(texImage[dev]);
        cudaCreateTextureObject(&texImage[dev], &texRes, &texDescr, NULL);
    }
}
//______________________________________________________________________________
//
//      Function:       matchedKernelParams
//
//      Description:    The 7 Point3D that kernelPixelBackprojection reads for projection i (geo: this piece of the
//                      image, with geo.alpha/theta/psi set for the projection). computeDeltasCube gives the voxel grid
//                      and the source in the scaled coordinate system rotated with the projection (detector rotation
//                      and offsets included); the rest is folded here, in double, so that the kernel's setup per
//                      projection is a few multiply-adds:
//                        [0..2] deltaX, deltaY, deltaZ: the change of the source-to-voxel vector per voxel index
//                        [3]    that vector at voxel (0,0,0): xyzOrigin - source, and COR (in detector pixels) in y
//                        [4]    (DSO-DSD-source.x, source.y+nDetecU/2, source.z+nDetecV/2): t=[4].x/vectX, u=vectY*t+[4].y
//                        [5]    (1/d_n, dDetecU^2, dDetecV^2) for the weight; d_n=source.x-(DSO-DSD) is the distance
//                               from the source to the detector plane
//                        [6]    unused
//______________________________________________________________________________

static void matchedKernelParams(const Geometry& geo, int i, Point3D* params)
{
    Point3D deltaX, deltaY, deltaZ, xyzOrigin, source;
    computeDeltasCube(geo, i, &xyzOrigin, &deltaX, &deltaY, &deltaZ, &source);
    const double COR=geo.COR[i], DSO=geo.DSO[i], DSD=geo.DSD[i];
    Point3D c0, tuv, wpar, unused={0,0,0};
    c0.x=(float)((double)xyzOrigin.x-(double)source.x);
    c0.y=(float)((double)xyzOrigin.y-COR/(double)geo.dDetecU-(double)source.y);
    c0.z=(float)((double)xyzOrigin.z-(double)source.z);
    tuv.x=(float)(DSO-DSD-(double)source.x);
    tuv.y=(float)((double)source.y+(double)geo.nDetecU*0.5);
    tuv.z=(float)((double)source.z+(double)geo.nDetecV*0.5);
    wpar.x=(float)(1.0/((double)source.x-(DSO-DSD)));
    wpar.y=(float)((double)geo.dDetecU*(double)geo.dDetecU);
    wpar.z=(float)((double)geo.dDetecV*(double)geo.dDetecV);
    params[0]=deltaX;
    params[1]=deltaY;
    params[2]=deltaZ;
    params[3]=c0;
    params[4]=tuv;
    params[5]=wpar;
    params[6]=unused;
}

#ifndef BACKPROJECTION_HPP
void splitCTbackprojection(const GpuIds& gpuids, Geometry geo,int nalpha, unsigned int* split_image, unsigned int * split_projections, bool* one_upload){
    
    
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
        // The smallest working set is one stack of projections and one image slice.
        const size_t mem_slice=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*sizeof(float);
        if (mem_GPU_global <= 2*mem_proj*PROJ_PER_KERNEL+mem_slice)
            tigreNotEnoughGpuMemory("Atb:voxel_backprojection2:GPUmemory","Backprojection",2*mem_proj*PROJ_PER_KERNEL+mem_slice,mem_GPU_global);
        size_t mem_free=mem_GPU_global-2*mem_proj*PROJ_PER_KERNEL;
        
        *split_image=(mem_image/deviceCount+mem_free-1)/mem_free;
        // Now knowing how many splits we have for images, we can recompute how many slices of projections actually
        // fit on the GPU. Must be more than 0 obviously.
        
        mem_free=mem_GPU_global-(mem_image/deviceCount)/(*split_image); // NOTE: There is some rounding error, but its in the order of bytes, and we have 5% of GPU free jsut in case. We are safe
        
        
        *split_projections=(mem_proj*PROJ_PER_KERNEL*2+mem_free-1)/mem_free;
        
    }
    // All the projections in one texture upload instead of chunks (thin detectors, see
    // tigre::projectionsInOneUpload): only with the whole problem in one piece, and if they fit next to the image.
    *one_upload = *split_image==1 && *split_projections==1 && tigre::projectionsInOneUpload(geo.nDetecV, nalpha) &&
                  mem_image/deviceCount+mem_proj*(size_t)nalpha < mem_GPU_global;
}


void computeDeltasCube(Geometry geo,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S)
{
    Point3Ddouble P, Px,Py,Pz;
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
    Point3Ddouble source;
    source.x=geo.DSD[i]; //already offseted for rotation
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
    
    
    *xyzorigin=P.to_float();
    *S=source.to_float();
}  // END computeDeltasCube

void checkFreeMemory(const GpuIds& gpuids,size_t *mem_GPU_global){
    size_t memfree;
    size_t memtotal;
    const int gpuids.GetLength();
    
    for (int dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(gpuids[dev]);
        cudaMemGetInfo(&memfree,&memtotal);
        if(dev==0) *mem_GPU_global=memfree;
        if(memfree<memtotal/2){
            mexWarnMsgIdAndTxt("voxel_backprojection:Atb:GPU","One (or more) of your GPUs is being heavily used by another program (less than half of its memory is free). TIGRE will try to run anyway, but it may need to split the problem into smaller pieces (slower) or run out of GPU memory.\n");
        }
        cudaCheckErrors("Check mem error");
        
        *mem_GPU_global=(memfree<*mem_GPU_global)?memfree:*mem_GPU_global;
    }
    *mem_GPU_global=(size_t)((double)*mem_GPU_global*0.95);
    
    //*mem_GPU_global= insert your known number here, in bytes.
}

#endif
