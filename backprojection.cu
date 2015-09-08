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
#include "projection.hpp"
#include "mex.h"
#include <math.h>

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            mexPrintf("%s \n",msg);\
            mexErrMsgIdAndTxt("CBCT:CUDA:interpolation",cudaGetErrorString(__err));\
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
    
___global___ void kernelPixelBackprojection(Geometry geo, double* image){
    
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ)
        return;
    
    int indVz = idx / (geo.nVoxelX*geo.nVoxelY);
    int resZ  = idx % (geo.nVoxelX*geo.nVoxelY);

    int indY= resZ / geo.nVoxelX;
    int indX= resZ % geo.nVoxelX;
    // Geometric trasnformations:
    
     Point3D S;
     S.x=geo.DSO;
     S.y=0;
     S.z=0;
     
     Point3D P;
     P.x=(indX-(geo.nVoxelX/2-0.5)) *geo.dVoxelX+geo.offOrigX;
     P.y=(indY-(geo.nVoxelY/2-0.5)) *geo.dVoxelY+geo.offOrigY;
     P.z=(indZ-(geo.nVoxelZ/2-0.5)) *geo.dVoxelZ+geo.offOrigZ;
     
     
             
     
    
    
    
    
    
    
    
    
    
    
}
    
    
int projection(float const * const projections, Geometry geo, double*** result,double const * const alphas,int nalpha){
 
    // BEFORE DOING ANYTHING: Use the proper CUDA enabled GPU: Tesla K40c
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
     if (deviceCount == 0)
    {
        mexErrMsgIdAndTxt("CBCT:CUDA:Ax:cudaGetDeviceCount","No CUDA enabled NVIDIA GPUs found");
    }
    bool found=false;
    for (int dev = 0; dev < deviceCount; ++dev)
    {
        cudaSetDevice(dev);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        if (strcmp(deviceProp.name, "Tesla K40c") == 0){
            found=true;
            break;
        }
    }
    if (!found)
        mexErrMsgIdAndTxt("CBCT:CUDA:Ax:cudaDevice","No Tesla K40c found");

    
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
    size_t num_bytes = geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ * sizeof(double);
    double* dProjection;
    cudaMalloc((void**)&dProjection, num_bytes);
    cudaMemset(dProjection,0,num_bytes);
    cudaCheckErrors("cudaMalloc fail");
    
    
}


