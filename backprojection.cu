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
#include "backprojection.hpp"
#include "mex.h"
#include <math.h>

// https://stackoverflow.com/questions/16282136/is-there-a-cuda-equivalent-of-perror
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            mexPrintf("%s \n",msg);\
            mexErrMsgIdAndTxt("CBCT:CUDA:backprojection",cudaGetErrorString(__err));\
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
    
__global__ void kernelPixelBackprojection(Geometry geo, 
                                            double* image,
                                            int indAlpha,
                                            Point3D deltaX ,
                                            Point3D deltaY, 
                                            Point3D deltaZ,
                                            Point3D xyzOrigin){
    //Make sure we dont go out of bounds
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nVoxelX* geo.nVoxelY *geo.nVoxelZ)
        return;
    
    int indZ = idx / (geo.nVoxelX*geo.nVoxelY);
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
     
     
     double t=(geo.DSO-geo.DSD /*-DDO*/ - S.x)/vectX;
     double y,z;
     y=vectY*t+S.y;
     z=vectZ*t+S.z;
     
     image[idx]+=tex3D(tex,(y-(geo.offDetecU/geo.dDetecU-(geo.nDetecU/2-0.5))) +0.5 ,
                           (z-(geo.offDetecV/geo.dDetecV-(geo.nDetecV/2-0.5))) +0.5 , 
                            indAlpha                                           +0.5);
     
     
//     image[idx]=-(geo.offDetecU/geo.dDetecU-(geo.nDetecU/2-0.5)) ;
    
}
    
    
int backprojection(float const * const projections, Geometry geo, double* result,double const * const alphas,int nalpha){
 
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
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        if (strcmp(deviceProp.name, "Tesla K40c") == 0){
            cudaSetDevice(dev);
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
        geo.alpha=alphas[i];
        computeDeltasCube(geo,geo.alpha,&xyzOrigin,&deltaX,&deltaY,&deltaZ);
        kernelPixelBackprojection<<<(geo.nVoxelX*geo.nVoxelY*geo.nVoxelZ + MAXTREADS-1) / MAXTREADS,MAXTREADS>>>
                (geo,dimage,i,deltaX,deltaY,deltaZ,xyzOrigin);
        cudaCheckErrors("Kernel fail");
    }
    cudaMemcpy(result, dimage, num_bytes, cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy fail");
    
     cudaUnbindTexture(tex);
     cudaCheckErrors("Unbind  fail");
     
     cudaFree(dimage);
     cudaFreeArray(d_projectiondata);
     cudaCheckErrors("cudaFree d_imagedata fail");
    return 0;
    
}
void computeDeltasCube(Geometry geo, double alpha, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ){
    
     Point3D P0, Px0,Py0,Pz0;
     // Get coords of Img(0,0,0)
     P0.x=-(geo.sVoxelX/2-geo.dVoxelX/2)-geo.offOrigX;
     P0.y=-(geo.sVoxelY/2-geo.dVoxelY/2)-geo.offOrigY;
     P0.z=-(geo.sVoxelZ/2-geo.dVoxelZ/2)-geo.offOrigZ;
     
     // Get coors from next voxel in each direction
     Px0.x=P0.x+geo.dVoxelX;       Py0.x=P0.x;                Pz0.x=P0.x;
     Px0.y=P0.y;                   Py0.y=P0.y+geo.dVoxelY;    Pz0.y=P0.y;
     Px0.z=P0.z;                   Py0.z=P0.z;                Pz0.z=P0.z+geo.dVoxelZ;
     
     // Rotate image (this is equivalent of rotating the source and detector)
     
     Point3D P, Px,Py,Pz; // We need other auxiliar variables to be able to perform the rotation, or we would overwrite values!
     P.x =P0.x *cos(alpha) -P0.y*sin(alpha);       P.y =P0.x *sin(alpha) +P0.y*cos(alpha);      P.z =P0.z;
     Px.x=Px0.x*cos(alpha)-Px0.y*sin(alpha);       Px.y=Px0.x*sin(alpha)+Px0.y*cos(alpha);      Px.z=Px0.z;
     Py.x=Py0.x*cos(alpha)-Py0.y*sin(alpha);       Py.y=Py0.x*sin(alpha)+Py0.y*cos(alpha);      Py.z=Py0.z;
     Pz.x=Pz0.x*cos(alpha)-Pz0.y*sin(alpha);       Pz.y=Pz0.x*sin(alpha)+Pz0.y*cos(alpha);      Pz.z=Pz0.z;
     
     // Scale coords so detector pixels are 1x1
     
     P.z =P.z /geo.dDetecV;          P.y =P.y /geo.dDetecU;
     Px.z=Px.z/geo.dDetecV;          Px.y=Px.y/geo.dDetecU;
     Py.z=Py.z/geo.dDetecV;          Py.y=Py.y/geo.dDetecU;
     Pz.z=Pz.z/geo.dDetecV;          Pz.y=Pz.y/geo.dDetecU;
     
     *xyzorigin=P;
     deltaX->x=Px.x-P.x;   deltaX->y=Px.y-P.y;    deltaX->z=Px.z-P.z;
     deltaY->x=Py.x-P.x;   deltaY->y=Py.y-P.y;    deltaY->z=Py.z-P.z;
     deltaZ->x=Pz.x-P.x;   deltaZ->y=Pz.y-P.y;    deltaZ->z=Pz.z-P.z;
     

}


