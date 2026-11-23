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
#include <vector>
#include "voxel_backprojection.hpp"
#include "voxel_backprojection_parallel.hpp"

#include "TIGRE_common.hpp"
#include "tigre_release_on_error.hpp"
#include "page_locking.hpp"
#include <math.h>

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
void CreateTextureParallel( float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream, bool allocate, cudaMemcpyKind kind);

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

// Per projection, the texture coordinates of the voxel (i,j,k) are linear in i, j and k (parallel rays: the voxel
// projects straight onto the detector): u = U.w + i*U.x + j*U.y + k*U.z and v = V.w + i*V.x + j*V.y + k*V.z,
// with U = projUVArrayDevParallel[2*p] and V = projUVArrayDevParallel[2*p+1], folded on the host (rotation,
// offsets, COR, detector rotation and the texel centre included).
__constant__ float4 projUVArrayDevParallel[2*PROJ_PER_KERNEL];

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
    unsigned long long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long long indX = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long startIndZ = blockIdx.z * VOXELS_PER_THREAD + threadIdx.z;  // the first voxel of this thread's column
    if (indX>=geo.nVoxelX || indY>=geo.nVoxelY || startIndZ>=geo.nVoxelZ)
        return;
    
    // The column of voxels this thread updates, kept in registers. Voxels of the column beyond the image (last block
    // in z) are computed too, harmlessly (a texture read), and never stored.
    float voxelColumn[VOXELS_PER_THREAD];
    const unsigned long long sliceSize=(unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY;
    const unsigned long long idx0=startIndZ*sliceSize+indY*(unsigned long long)geo.nVoxelX+indX;
#pragma unroll
    for(int colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
        voxelColumn[colIdx] = (startIndZ+colIdx<geo.nVoxelZ) ? image[idx0+colIdx*sliceSize] : 0.0f;
    
    const float fIndX=(float)indX, fIndY=(float)indY, fStartIndZ=(float)startIndZ;
    for(int projNumber=0; projNumber<PROJ_PER_KERNEL; projNumber++)
    {
        const int indAlpha = currProjSetNumber*PROJ_PER_KERNEL+projNumber;  // the projection's index in the texture
        if(indAlpha>=totalNoOfProjections)
            break;
        const float4 U = projUVArrayDevParallel[2*projNumber];
        const float4 V = projUVArrayDevParallel[2*projNumber+1];
        // The first voxel of the column, then one step of (U.z, V.z) per voxel. Before, every voxel rebuilt its
        // position from the indices (9 multiply-adds), the source on its ray and the ray's intersection with the
        // detector (two divisions): for parallel rays these cancel, the detector coordinates are the voxel's own.
        float u=U.w+fIndX*U.x+fIndY*U.y+fStartIndZ*U.z;
        float v=V.w+fIndX*V.x+fIndY*V.y+fStartIndZ*V.z;
        const int layer=(int)indAlpha;  // the projection: a layer of the layered texture
#pragma unroll
        for(int colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
        {
#if IS_FOR_MATLAB_TIGRE
            voxelColumn[colIdx]+=tex2DLayered<float>(tex, v, u, layer);
#else
            voxelColumn[colIdx]+=tex2DLayered<float>(tex, u, v, layer);
#endif
            u+=U.z; v+=V.z;
        }
    }
    
#pragma unroll
    for(int colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
        if (startIndZ+colIdx<geo.nVoxelZ)
            image[idx0+colIdx*sliceSize] = voxelColumn[colIdx];  // = and not +=: the column was read at the start
    
}  // END kernelPixelBackprojectionFDK




//______________________________________________________________________________
//
//      Function:       voxel_backprojection_parallel
//
//      Description:    Main host function for FDK backprojection (invokes the kernel)
//______________________________________________________________________________

int voxel_backprojection_parallel(float  *  projections, Geometry geo, float* result,float const * const alphas, int nalpha, const GpuIds& gpuids, bool gpuData)
{
    // The image is split in z between the GPUs (and, if a GPU cannot hold its part, in more pieces, done one after
    // the other). Each GPU backprojects all the projections into its piece, PROJ_PER_KERNEL projections at a time
    // through two texture arrays, so that the next block's upload overlaps this block's kernel (#286).
    const int deviceCount = gpuids.GetLength();
    // If a CUDA error raises (Python build) before the end of this function,
    // onError releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Atb:voxel_backprojection_parallel:GPUselect","There are no available device(s) that support CUDA\n");
    }
    if (gpuData && deviceCount != 1) {
        mexErrMsgIdAndTxt("Atb:GPUdata", "GPU arrays can only be backprojected on one GPU, the one they are on.");
    }
    const cudaMemcpyKind kindIn  = gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice;
    const cudaMemcpyKind kindOut = gpuData ? cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost;
    std::vector<int> devIds(deviceCount);  // copied: the release actions must not hold on to gpuids
    for (int d = 0; d < deviceCount; d++) { devIds[d] = gpuids[d]; onError.addDevice(devIds[d]); }
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Atb:voxel_backprojection_parallel:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;

    const unsigned int nBlocks=(nalpha+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;  // blocks of PROJ_PER_KERNEL projections
    const int nArrays=(nBlocks>1)? 2 : 1;
#if IS_FOR_MATLAB_TIGRE
    const cudaExtent extentProj=make_cudaExtent(geo.nDetecV, geo.nDetecU, (nalpha<PROJ_PER_KERNEL)? nalpha : PROJ_PER_KERNEL);
#else
    const cudaExtent extentProj=make_cudaExtent(geo.nDetecU, geo.nDetecV, (nalpha<PROJ_PER_KERNEL)? nalpha : PROJ_PER_KERNEL);
#endif
    const size_t mem_slice=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*sizeof(float);
    const size_t mem_arrays=(size_t)nArrays*extentProj.width*extentProj.height*extentProj.depth*sizeof(float);

    // How many pieces per GPU? Estimate from the free memory, then allocate: the image piece and the projection
    // arrays on every GPU. If that runs out of memory (3D arrays need more than their raw size), split more and try
    // again, down to one slice per piece.
    size_t mem_GPU_global=0;
    for (dev = 0; dev < deviceCount; dev++){
        size_t memfree=0, memtotal=0;
        cudaSetDevice(devIds[dev]);
        cudaMemGetInfo(&memfree,&memtotal);
        if (dev==0 || memfree<mem_GPU_global) mem_GPU_global=memfree;
    }
    mem_GPU_global=(size_t)((double)mem_GPU_global*0.95);
    cudaCheckErrors("Check mem error");
    const unsigned int maxPieces=(unsigned int)geo.nVoxelZ;  // pieces in total, over all GPUs
    unsigned int splits=1;  // pieces per GPU
    if (mem_GPU_global <= mem_arrays+mem_slice)
        tigreNotEnoughGpuMemory("Atb:voxel_backprojection_parallel:GPUmemory","Backprojection (parallel beam)",mem_arrays+mem_slice,mem_GPU_global);
    {
        const size_t mem_piece=mem_slice*(size_t)((geo.nVoxelZ+deviceCount-1)/deviceCount);
        const size_t mem_free=mem_GPU_global-mem_arrays;
        if (mem_piece>=mem_free) splits=(unsigned int)(mem_piece/mem_free+1);
        while (splits>1 && (unsigned long long)splits*deviceCount>maxPieces) splits--;
    }
    // GPU data, whole image in one piece: backproject straight into the caller's array.
    float** dimage=(float**)malloc(deviceCount*sizeof(float*));
    cudaArray** d_cuArrTex=(cudaArray**)malloc(2*deviceCount*sizeof(cudaArray*));
    Geometry* geoArray=NULL;
    bool directOut=false;
    while (true){
        memset(dimage, 0, deviceCount*sizeof(float*));
        memset(d_cuArrTex, 0, 2*deviceCount*sizeof(cudaArray*));
        geoArray=(Geometry*)malloc((size_t)splits*deviceCount*sizeof(Geometry));
        createGeoArray(splits*deviceCount,geo,geoArray,nalpha);
        directOut = gpuData && splits==1 && deviceCount==1;
        const size_t num_bytes_img=mem_slice*(size_t)geoArray[0].nVoxelZ;
        cudaError_t err=cudaSuccess;
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        for (dev = 0; dev < deviceCount && err==cudaSuccess; dev++){
            cudaSetDevice(devIds[dev]);
            if (directOut) dimage[dev]=result;
            else           err=cudaMalloc((void**)&dimage[dev], num_bytes_img);
            for (int a = 0; a < nArrays && err==cudaSuccess; a++)
                err=cudaMalloc3DArray(&d_cuArrTex[2*dev+a], &channelDesc, extentProj, cudaArrayLayered);
        }
        if (err==cudaSuccess) break;
        cudaGetLastError();
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(devIds[dev]);
            if (dimage[dev] && dimage[dev]!=result) cudaFree(dimage[dev]);
            for (int a = 0; a < 2; a++) if (d_cuArrTex[2*dev+a]) cudaFreeArray(d_cuArrTex[2*dev+a]);
        }
        cudaGetLastError();
        const unsigned int piece=geoArray[0].nVoxelZ;
        freeGeoArray(splits*deviceCount, geoArray);
        if (err!=cudaErrorMemoryAllocation){
            free(dimage); free(d_cuArrTex);
            tigreCudaError("Atb:voxel_backprojection_parallel:GPUmemory", "Malloc error", cudaGetErrorString(err));
        }
        if (piece<=1 || (unsigned long long)splits*deviceCount>=maxPieces){
            free(dimage); free(d_cuArrTex);
            size_t free_bytes=0, total_bytes=0;
            cudaMemGetInfo(&free_bytes, &total_bytes);
            tigreNotEnoughGpuMemory("Atb:voxel_backprojection_parallel:GPUmemory","Backprojection (parallel beam)",mem_arrays+mem_slice,free_bytes);
        }
        // smaller pieces
        while ((unsigned long long)splits*deviceCount<maxPieces && (geo.nVoxelZ+splits*deviceCount-1)/(splits*deviceCount)>=piece) splits++;
    }
    const unsigned int nPieces=splits*deviceCount;
    onError.add([geoArray, nPieces]{ freeGeoArray(nPieces, geoArray); });
    onError.add([dimage, d_cuArrTex, deviceCount, devIds, result]{
        for (int d = 0; d < deviceCount; d++) {
            cudaSetDevice(devIds[d]);
            if (dimage[d] && dimage[d] != result) cudaFree(dimage[d]);
            for (int a = 0; a < 2; a++) if (d_cuArrTex[2*d+a]) cudaFreeArray(d_cuArrTex[2*d+a]);
        }
        free(dimage); free(d_cuArrTex);
    });

    // Two streams per GPU: kernels on the first, projection uploads on the second.
    const int nStreams=2*deviceCount;
    cudaStream_t* stream=(cudaStream_t*)malloc(nStreams*sizeof(cudaStream_t));
    memset(stream, 0, nStreams*sizeof(cudaStream_t));
    onError.add([stream, nStreams]{ for (int s = 0; s < nStreams; s++) if (stream[s]) cudaStreamDestroy(stream[s]); free(stream); });
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(devIds[dev]);
        for (int i = 0; i < 2; ++i)
            cudaStreamCreate(&stream[i+dev*2]);
    }
    cudaCheckErrors("Stream creation fail");

    // Page-lock the projections and the result when the image is split or several GPUs are used (tigre::worthPageLocking):
    // then the uploads are sent again for each piece or to each GPU, and the GPUs copy their pieces of the result back at
    // the same time instead of one after another.
    bool registered=false, result_registered=false;
    int isHostRegisterSupported = 0;
#if CUDART_VERSION >= 9020
    cudaDeviceGetAttribute(&isHostRegisterSupported,cudaDevAttrHostRegisterSupported,devIds[0]);
#endif
    if (tigre::worthPageLocking(isHostRegisterSupported, gpuData, splits, deviceCount)){
        registered=tigre::hostRegister(projections, (size_t)geo.nDetecU*(size_t)geo.nDetecV*(size_t)nalpha*sizeof(float));
        // The caller's array: left registered, the next call with it would fail.
        if (registered) onError.add([projections]{ cudaHostUnregister(projections); });
        result_registered=tigre::hostRegister(result, (size_t)geo.nVoxelX*(size_t)geo.nVoxelY*(size_t)geo.nVoxelZ*sizeof(float));
        if (result_registered) onError.add([result]{ cudaHostUnregister(result); });
    }
    cudaCheckErrors("Error pinning memory");
    // The geometry constants of each launch are filled here and copied to __constant__ memory before the kernel.
    // Ordinary host memory: page-locking these few KB (allocating and freeing them with cudaMallocHost and
    // cudaFreeHost) cost ~0.45 ms per call and did not change the kernels' time (measured on an RTX 4070: the
    // kernels read the __constant__ copy either way).

    float4* projUVArrayHostParallel=(float4*)malloc(2*PROJ_PER_KERNEL*sizeof(float4));
    if (!projUVArrayHostParallel) mexErrMsgIdAndTxt("TIGRE:hostMemory","Not enough host memory for the geometry constants");
    onError.add([projUVArrayHostParallel]{ free(projUVArrayHostParallel); });

    // Texture objects, two per GPU. Zeroed so CreateTextureParallel can tell an empty slot from one holding a
    // texture object it must destroy first.
    cudaTextureObject_t *texProj = new cudaTextureObject_t[2*deviceCount]();
    onError.add([texProj, deviceCount, devIds]{
        for (int s = 0; s < 2*deviceCount; s++) {
            cudaSetDevice(devIds[s/2]);
            if (texProj[s]) cudaDestroyTextureObject(texProj[s]);
        }
        delete[] texProj;
    });

    const int divx=16, divy=32, divz=VOXELS_PER_THREAD;  // see Zinsser and Keck; each thread does a column of divz voxels
    for (unsigned int sp=0; sp<splits; sp++){
        // GPU dev works on piece sp*deviceCount+dev (some trailing pieces may be empty).
        for (dev = 0; dev < deviceCount; dev++){
            const Geometry& g=geoArray[sp*deviceCount+dev];
            if (g.nVoxelZ<=0) continue;
            cudaSetDevice(devIds[dev]);
            cudaMemsetAsync(dimage[dev],0,mem_slice*(size_t)g.nVoxelZ,stream[dev*2]);
        }
        for (unsigned int b=0; b<nBlocks; b++){
            const unsigned int first=b*PROJ_PER_KERNEL;
            const unsigned int nproj=((unsigned int)nalpha-first<PROJ_PER_KERNEL)? (unsigned int)nalpha-first : PROJ_PER_KERNEL;
            // Upload this block of projections to every GPU working on a piece.
            for (dev = 0; dev < deviceCount; dev++){
                if (geoArray[sp*deviceCount+dev].nVoxelZ<=0) continue;
                cudaSetDevice(devIds[dev]);
                CreateTextureParallel(&projections[(size_t)first*geo.nDetecU*geo.nDetecV],geo,&d_cuArrTex[2*dev+(b%2)],nproj,&texProj[2*dev+(b%2)],&stream[dev*2],false,kindIn);
            }
            for (dev = 0; dev < deviceCount; dev++){
                Geometry& g=geoArray[sp*deviceCount+dev];
                if (g.nVoxelZ<=0) continue;
                cudaSetDevice(devIds[dev]);
                for (unsigned int j=0; j<nproj; j++){
                    const unsigned int p=first+j;  // global projection index
                    Point3D deltaX,deltaY,deltaZ,xyzOrigin,source;
                    g.alpha=-alphas[p*3];
                    g.theta=-alphas[p*3+1];
                    g.psi  =-alphas[p*3+2];
                    computeDeltasCubeParallel(g,p,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
                    // The kernel's old per-voxel sum, folded: u = P.y - COR/dDetecU + nDetecU/2 - 0.5 and
                    // v = P.z + nDetecV/2 - 0.5 for the voxel's position P = xyzOrigin + i*deltaX + j*deltaY +
                    // k*deltaZ (detector units), plus 0.5 for the texel centre.
                    const double u0=(double)xyzOrigin.y-(double)g.COR[p]/(double)geo.dDetecU+geo.nDetecU/2.0;
                    const double v0=(double)xyzOrigin.z+geo.nDetecV/2.0;
                    projUVArrayHostParallel[2*j]  =make_float4(deltaX.y,deltaY.y,deltaZ.y,(float)u0);
                    projUVArrayHostParallel[2*j+1]=make_float4(deltaX.z,deltaY.z,deltaZ.z,(float)v0);
                }
                cudaMemcpyToSymbolAsync(projUVArrayDevParallel, projUVArrayHostParallel, sizeof(float4)*2*PROJ_PER_KERNEL,0,cudaMemcpyHostToDevice,stream[dev*2]);
                // The host arrays are refilled for the next GPU; and this also waits for the kernel of the block
                // before, which read the texture array the next block's upload overwrites.
                cudaStreamSynchronize(stream[dev*2]);
                cudaStreamSynchronize(stream[dev*2+1]);  // this block's upload
                dim3 grid((g.nVoxelX+divx-1)/divx, (g.nVoxelY+divy-1)/divy, (g.nVoxelZ+divz-1)/divz);
                dim3 block(divx,divy,1);  // 1 in z: each thread works on a column of VOXELS_PER_THREAD voxels
                kernelPixelBackprojection_parallel<<<grid,block,0,stream[dev*2]>>>(g,dimage[dev],0,nproj,texProj[2*dev+(b%2)]);
            }
            cudaCheckErrors("Kernel fail");
        }
        // Copy the pieces out.
        for (dev = 0; dev < deviceCount; dev++){
            const unsigned int p=sp*deviceCount+dev;
            if (geoArray[p].nVoxelZ<=0) continue;
            cudaSetDevice(devIds[dev]);
            if (!directOut)
                cudaMemcpyAsync(&result[mem_slice/sizeof(float)*(size_t)p*(size_t)geoArray[0].nVoxelZ], dimage[dev], mem_slice*(size_t)geoArray[p].nVoxelZ, kindOut, stream[dev*2]);
        }
        for (dev = 0; dev < deviceCount; dev++){
            cudaSetDevice(devIds[dev]);
            cudaDeviceSynchronize();
        }
        cudaCheckErrors("cudaMemcpy result fail");
    }

    onError.dismiss();  // from here on the function releases everything itself
    for (dev = 0; dev < deviceCount; dev++){
        cudaSetDevice(devIds[dev]);
        for (int a = 0; a < 2; a++){
            if (texProj[2*dev+a]) cudaDestroyTextureObject(texProj[2*dev+a]);
            if (d_cuArrTex[2*dev+a]) cudaFreeArray(d_cuArrTex[2*dev+a]);
        }
        if (!directOut) cudaFree(dimage[dev]);
        cudaStreamDestroy(stream[dev*2]);
        cudaStreamDestroy(stream[dev*2+1]);
    }
    delete[] texProj;
    free(d_cuArrTex);
    free(dimage);
    free(stream);
    freeGeoArray(nPieces, geoArray);
    free(projUVArrayHostParallel);
    if (registered) cudaHostUnregister(projections);
    if (result_registered) cudaHostUnregister(result);
    cudaCheckErrors("cudaFree fail");
    return 0;
}  // END voxel_backprojection_parallel


void computeDeltasCubeParallel(Geometry geo, int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D *S)
{
    
    Point3Ddouble P, Px,Py,Pz;
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
    
    
    Point3Ddouble source;
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
    
    
    // cast the results from the double precision calculations back to float
    *xyzorigin=P.to_float();
    *S=source.to_float();

    
}  // END computeDeltasCube
void CreateTextureParallel(float* projectiondata,Geometry geo,cudaArray** d_cuArrTex,unsigned int nangles, cudaTextureObject_t *texImage,cudaStream_t* stream, bool alloc, cudaMemcpyKind kind)
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
        tigre::checkMalloc(cudaMalloc3DArray(&d_cuArrTex[0], &channelDesc, extent, cudaArrayLayered), "Atb:voxel_backprojection_parallel:GPUmemory", "Backprojection (parallel beam)", (size_t)extent.width*extent.height*extent.depth*sizeof(float));
        cudaCheckErrors("Texture memory allocation fail");
        }
        cudaMemcpy3DParms copyParams = {0};
        
        
        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)projectiondata, extent.width*sizeof(float), extent.width, extent.height);
        copyParams.dstArray = d_cuArrTex[0];
        copyParams.extent   = extent;
        copyParams.kind     = kind;
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
        // Each projection block (or image split) creates a new texture object in this
        // slot. Destroy the one it replaces, otherwise every call leaks one per block
        // until cudaCreateTextureObject fails with "out of memory" (#617). The kernels
        // that used it have finished: its array has just been overwritten above.
        if (texImage[0]) cudaDestroyTextureObject(texImage[0]);
        cudaCreateTextureObject(&texImage[0], &texRes, &texDescr, NULL);
        cudaCheckErrors("Texture object creation fail");
    
}