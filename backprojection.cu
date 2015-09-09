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
    
___global___ void kernelPixelBackprojection(Geometry geo, 
                                            double* image,
                                            int indAlpha,
                                            Point3D deltaX ,
                                            Point3D deltaY, 
                                            Point3D deltaZ,
                                            Point3D xyzOrigin){){
    //Make sure we dont go out of bounds
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ)
        return;
    
    int indz = idx / (geo.nVoxelX*geo.nVoxelY);
    int resZ  = idx % (geo.nVoxelX*geo.nVoxelY);

    int indY= resZ / geo.nVoxelX;
    int indX= resZ % geo.nVoxelX;
    // Geometric trasnformations:
    
    //Source
     Point3D S;   
     S.x=geo.DSO;
     S.y=geo.offDetecU/geo.dDetecU;
     S.z=geo.offDetecV/geo.dDetecV;
     
     Point3D P;
     P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
     P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y);
     P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
     
     double vectX,vectY,vectZ;
     vectX=(P.x -S.x); 
     vectY=(P.y -S.y); 
     vectZ=(P.z -S.z); 
     
     
     double t=(geo.DSD-geo.DSO /*DDO*/ - S.x)/vectx;
     double y,z;
     y=vectY*t+S.y;
     z=vectZ*t+S.z;
     
     image[idx]+=tex3D(tex,(y+(geo.offDetecU-(geo.sDetecU/2-0.5)))/geo.dDetecU +0.5 ,
                           (z+(geo.offDetecV-(geo.sDetecV/2-0.5)))/geo.dDetecV +0.5 , 
                            indAlpha                                           +0.5);
     
    
    
    
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
    double* dimage;
    cudaMalloc((void**)&dimage, num_bytes);
    cudaMemset(dimage,0,num_bytes);
    cudaCheckErrors("cudaMalloc fail");
    
    
    
    Point3D deltaX,deltaY,deltaZ,xyzOrigin;
    
    for (int i=0;i<nalpha;i++){
        
        
        kernelPixelBackprojection<<<(geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ + MAXTREADS-1) / MAXTREADS,MAXTREADS>>>
                (geo,dimage,i,deltaX,deltaY,deltaZ,xyzOrigin);
        
    }
    
}
void computeDeltasCube(Geometry geo, double alpha, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ){


