/*-------------------------------------------------------------------------
 *
 * SYCL functions for ray-voxel intersection based projection
 *
 * This file has the necessary fucntiosn to perform X-ray parallel projection
 * operation given a geaometry, angles and image. It usesthe so-called
 * Jacobs algorithm to compute efficiently the length of the x-rays over
 * voxel space. Its called Siddon because Jacobs algorithm its just a small
 * improvement over the traditional Siddons method.
 *
 * CODE by       Ander Biguri
 *
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
#include <vector>
#include <sycl/sycl.hpp>
#include "Siddon_projection.hpp"
#include "Siddon_projection_parallel.hpp"
#include "TIGRE_common.hpp"
#include "tigre_sycl.hpp"
#include "tigre_release_on_error.hpp"
#include <math.h>

#define MAXTREADS 1024
#define PROJ_PER_BLOCK 9
#define PIXEL_SIZE_BLOCK 9
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
 *    --->y   |   |        | /              |
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



// Parameters of the PROJ_PER_BLOCK projections of one kernel call (constant memory in the CUDA code)
struct ProjParamsParallel { Point3D data[4*PROJ_PER_BLOCK]; };

static inline void kernelPixelDetector_parallel(const Geometry& geo,
        float* detector, const int currProjSetNumber, const int totalNoOfProjections, const float* tex,
        const Point3D* projParamsArrayDev, sycl::id<3> id){
    
    unsigned long long u = id[2];
    unsigned long long v = id[1];
    unsigned long long projNumber=id[0];
            
    if (u>= geo.nDetecU || v>= geo.nDetecV || projNumber>=PROJ_PER_BLOCK)
        return;
    
    unsigned long indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    
#if IS_FOR_MATLAB_TIGRE
    size_t idx =  (size_t)(u  * (unsigned long long)geo.nDetecV + v)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#else
    size_t idx =  (size_t)(v  * (unsigned long long)geo.nDetecU + u)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#endif
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    

    /////// Get coordinates XYZ of pixel UV
    unsigned long pixelV = geo.nDetecV-v-1;
    unsigned long pixelU = u;
    Point3D pixel1D;
    pixel1D.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    pixel1D.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    pixel1D.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
    
    
    source.x=(source.x+pixelU*deltaU.x+pixelV*deltaV.x);
    source.y=(source.y+pixelU*deltaU.y+pixelV*deltaV.y);
    source.z=(source.z+pixelU*deltaU.z+pixelV*deltaV.z);
    ///////
    // Siddon's ray-voxel intersection, optimized as in doi=10.1.1.55.7516
    //////
    Point3D ray;
    // vector of Xray
    ray.x=pixel1D.x-source.x;
    ray.y=pixel1D.y-source.y;
    ray.z=pixel1D.z-source.z;

    float axm, aym, azm;
    float axM, ayM, azM;

    if (ray.x == 0.0f) {
        if (source.x >= 0.0f && source.x < geo.nVoxelX) {
            axm = -INFINITY;
            axM = INFINITY;
        } else {
            axm = INFINITY;
            axM = -INFINITY;
        }
    } else {
        float t1 = -source.x / ray.x;
        float t2 = (geo.nVoxelX - source.x) / ray.x;
        axm = sycl::fmin(t1, t2);
        axM = sycl::fmax(t1, t2);
    }

    if (ray.y == 0.0f) {
        if (source.y >= 0.0f && source.y < geo.nVoxelY) {
            aym = -INFINITY;
            ayM = INFINITY;
        } else {
            aym = INFINITY;
            ayM = -INFINITY;
        }
    } else {
        float t1 = -source.y / ray.y;
        float t2 = (geo.nVoxelY - source.y) / ray.y;
        aym = sycl::fmin(t1, t2);
        ayM = sycl::fmax(t1, t2);
    }

    if (ray.z == 0.0f) {
        if (source.z >= 0.0f && source.z < geo.nVoxelZ) {
            azm = -INFINITY;
            azM = INFINITY;
        } else {
            azm = INFINITY;
            azM = -INFINITY;
        }
    } else {
        float t1 = -source.z / ray.z;
        float t2 = (geo.nVoxelZ - source.z) / ray.z;
        azm = sycl::fmin(t1, t2);
        azM = sycl::fmax(t1, t2);
    }

    float am = sycl::fmax(sycl::fmax(axm, aym), azm);
    float aM = sycl::fmin(sycl::fmin(axM, ayM), azM);
    
    // line intersects voxel space ->   am<aM
    // See Siddon_projection.cu: without the return a missed ray falls
    // through to an Np built from indices that only mean something for a
    // hit; for an axis-parallel ray with the source outside the volume that
    // is (-inf)*0 = NaN, and (unsigned long)NaN = 2^63 spins forever. The
    // hang was NOT reproduced in this kernel (26 configurations tried); the
    // return here is by inspection.
    if (am>=aM){
        detector[idx]=0.0f;
        return;
    }

    // Walk the voxels the ray crosses, from where it enters the box (am) to where
    // it leaves (aM), adding up length x value. Same traversal as the cone-beam
    // kernel (Siddon_projection.cpp, #859): the first voxel is the one the ray is
    // in just after entering, the next crossing on each axis is that voxel's own
    // face ahead, and each crossing is computed from its plane.
    float pX=source.x+am*ray.x;
    float pY=source.y+am*ray.y;
    float pZ=source.z+am*ray.z;
    int i=(ray.x<0.0f)? (int)sycl::ceil(pX)-1 : (int)sycl::floor(pX);
    int j=(ray.y<0.0f)? (int)sycl::ceil(pY)-1 : (int)sycl::floor(pY);
    int k=(ray.z<0.0f)? (int)sycl::ceil(pZ)-1 : (int)sycl::floor(pZ);
    i=(i<0)? 0 : ((i>(int)geo.nVoxelX-1)? (int)geo.nVoxelX-1 : i);
    j=(j<0)? 0 : ((j>(int)geo.nVoxelY-1)? (int)geo.nVoxelY-1 : j);
    k=(k<0)? 0 : ((k>(int)geo.nVoxelZ-1)? (int)geo.nVoxelZ-1 : k);

    // Direction of update (eq 29), and the coordinate of the face ahead on each
    // axis: index+1 moving up, index moving down. The voxel is kept as integer
    // indices, which the fetch takes; the face coordinates are kept as floats
    // (exact: integers) so no conversion is needed in the loop.
    int sI=(ray.x<0.0f)? -1 : 1;
    int sJ=(ray.y<0.0f)? -1 : 1;
    int sK=(ray.z<0.0f)? -1 : 1;
    float fX=(float)(i+(sI>0));
    float fY=(float)(j+(sJ>0));
    float fZ=(float)(k+(sK>0));

    // Parameter of the next voxel boundary on each axis (infinite if the ray is
    // parallel to it): the boundary at coordinate b is at b*rX+oX. Each one is
    // computed from its plane, not by adding a step to the previous one (eq 28),
    // which drifts in float32 over the many boundaries of a large image.
    float rX=(ray.x!=0.0f)? 1.0f/ray.x : 0.0f;
    float rY=(ray.y!=0.0f)? 1.0f/ray.y : 0.0f;
    float rZ=(ray.z!=0.0f)? 1.0f/ray.z : 0.0f;
    float oX=-source.x*rX;
    float oY=-source.y*rY;
    float oZ=-source.z*rZ;
    float ax=(ray.x!=0.0f)? sycl::fma(fX,rX,oX) : INFINITY;
    float ay=(ray.y!=0.0f)? sycl::fma(fY,rY,oY) : INFINITY;
    float az=(ray.z!=0.0f)? sycl::fma(fZ,rZ,oZ) : INFINITY;

    float maxlength = sycl::sqrt(ray.x * ray.x * geo.dVoxelX * geo.dVoxelX + ray.y * ray.y * geo.dVoxelY * geo.dVoxelY + ray.z * ray.z * geo.dVoxelZ * geo.dVoxelZ);
    float sum=0.0f;
    float ac=am;
    // At most one step per voxel boundary, plus slack: the loop always ends, even
    // if rounding never quite lets a crossing reach aM. One fetch per step, before
    // the branch on the axis crossed.
    const int maxsteps=(int)geo.nVoxelX+(int)geo.nVoxelY+(int)geo.nVoxelZ+3;
    for (int step=0; step<maxsteps; step++){
        float aminc=sycl::fmin(sycl::fmin(ax,ay),az);
        sum+=(sycl::fmin(aminc,aM)-ac)*tigreFetchPoint(tex, geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ, i, j, k);
        if (aminc>=aM)
            break;
        if (ax==aminc){
            i+=sI; fX+=(float)sI; ac=ax; ax=sycl::fma(fX,rX,oX);
        }else if(ay==aminc){
            j+=sJ; fY+=(float)sJ; ac=ay; ay=sycl::fma(fY,rY,oY);
        }else{
            k+=sK; fZ+=(float)sK; ac=az; az=sycl::fma(fZ,rZ,oZ);
        }
    }
    detector[idx]=maxlength*sum;
}


int siddon_ray_projection_parallel(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // GPU arrays (#857): the pointers would be CUDA device memory, which a SYCL device
    // cannot read. tigre.Ax and tigre.Atb pass such arrays through the host instead.
    if (gpuData) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection_parallel:GPUdata","GPU arrays (PyTorch CUDA tensors, CuPy arrays) are not supported by the SYCL build: pass host (NumPy) arrays.\n");
    }
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection_parallel:GPUselect","There are no available SYCL device(s)\n");
    }
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Ax:Siddon_projection_parallel:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;
    // If an error raises (Python build) before the end of this function, onError
    // releases what it holds; see tigre_release_on_error.hpp.
    tigre::ReleaseOnError onError;
    try {
        // Check the available device memory: same image split rule as the CUDA code.
        size_t mem_GPU_global;
        checkFreeMemory(gpuids, &mem_GPU_global);

        size_t mem_image=                 (unsigned long long)geo.nVoxelX*(unsigned long long)geo.nVoxelY*(unsigned long long)geo.nVoxelZ*sizeof(float);
        size_t mem_proj=                  (unsigned long long)geo.nDetecU*(unsigned long long)geo.nDetecV*sizeof(float);

        // Does everything fit in the device?
        const bool fits_in_memory = mem_image+2*PROJ_PER_BLOCK*mem_proj<mem_GPU_global;
        unsigned int splits=1;
        if (!fits_in_memory) {
            // Nope nope.
            // approx free memory we have. We already have left some extra 5% free for internal stuff
            // we need a second projection memory to combine multi-GPU stuff.
            // The smallest working set is the 4 projection buffers and one image slice.
            const size_t mem_slice=(size_t)geo.nVoxelX*(size_t)geo.nVoxelY*sizeof(float);
            if (mem_GPU_global <= 4*PROJ_PER_BLOCK*mem_proj+mem_slice)
                tigreNotEnoughGpuMemory("Ax:Siddon_projection_parallel:GPUmemory","Forward projection (Siddon, parallel beam)",4*PROJ_PER_BLOCK*mem_proj+mem_slice,mem_GPU_global);
            size_t mem_free=mem_GPU_global-4*PROJ_PER_BLOCK*mem_proj;
            splits=mem_image/mem_free+1;// Ceil of the truncation
            if (splits>geo.nVoxelZ) splits=geo.nVoxelZ;  // at least one slice per piece
            // The pieces have ceil(nVoxelZ/splits) slices. Fewer pieces of that size can
            // cover the image (e.g. 64 slices in 12 pieces: 11 of 6 slices); the pieces
            // beyond would have a negative (wrapped around) number of slices.
            const unsigned int slices_per_piece=(geo.nVoxelZ+splits-1)/splits;
            splits=(geo.nVoxelZ+slices_per_piece-1)/slices_per_piece;
        }
        Geometry* geoArray = (Geometry*)malloc(splits*sizeof(Geometry));
        splitImage(splits,geo,geoArray,nangles);
        onError.add([geoArray, splits]{ freeGeoArray(splits, geoArray); });

        // Allocate the image (a slab of it if split) and one block of projections per device
        size_t num_voxels_split = (size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        size_t num_pixels_block = (size_t)PROJ_PER_BLOCK*(size_t)geo.nDetecU*(size_t)geo.nDetecV;
        std::vector<float*> dImage(deviceCount), dProjection(deviceCount);
        for (dev = 0; dev < deviceCount; dev++){
            sycl::queue* q = &tigreSyclQueue(gpuids[dev]);
            onError.addQueue(q);
            float* image_d = tigre::deviceMalloc<float>(num_voxels_split, gpuids[dev], "Ax:Siddon_projection_parallel:malloc", "Forward projection (Siddon, parallel beam)");
            onError.add([image_d, q]{ sycl::free(image_d, *q); });
            float* projection_d = tigre::deviceMalloc<float>(num_pixels_block, gpuids[dev], "Ax:Siddon_projection_parallel:malloc", "Forward projection (Siddon, parallel beam)");
            onError.add([projection_d, q]{ sycl::free(projection_d, *q); });
            dImage[dev] = image_d;
            dProjection[dev] = projection_d;
        }
        std::vector<float> partial(splits>1 ? num_pixels_block*deviceCount : 0);

        // Split the angles between devices
        int nangles_device=(nangles+deviceCount-1)/deviceCount;
        unsigned int noOfKernelCalls = (nangles_device+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_BLOCK

        ProjParamsParallel params;
        for (unsigned int sp=0;sp<splits;sp++){
            // Load the image (or the slab of it) on every device
            size_t linear_idx_start= (size_t)sp*(size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
            size_t num_voxels_curr = (size_t)geoArray[sp].nVoxelX*(size_t)geoArray[sp].nVoxelY*(size_t)geoArray[sp].nVoxelZ;
            for (dev = 0; dev < deviceCount; dev++)
                tigreSyclQueue(gpuids[dev]).memcpy(dImage[dev], &img[linear_idx_start], num_voxels_curr*sizeof(float));

            const Geometry geoSplit = geoArray[sp];
            for (unsigned int i=0; i<noOfKernelCalls; i++) {
                std::vector<sycl::event> copies(deviceCount);
                std::vector<int> projection_this_block(deviceCount, 0);
                for (dev=0;dev<deviceCount;dev++){
                    unsigned int proj_start=i*PROJ_PER_BLOCK+dev*nangles_device;
                    if (proj_start>=nangles || i*PROJ_PER_BLOCK>=nangles_device)
                        continue;
                    projection_this_block[dev]=std::min(std::min(PROJ_PER_BLOCK, nangles_device-(int)(i*PROJ_PER_BLOCK)), nangles-(int)proj_start);
                    for(unsigned int j=0; j<PROJ_PER_BLOCK; j++){
                        unsigned int proj_global=(i*PROJ_PER_BLOCK+j)+dev*nangles_device;
                        if (proj_global>=nangles)
                            break;
                        if ((i*PROJ_PER_BLOCK+j)>=nangles_device)
                            break;
                        geoArray[sp].alpha=angles[proj_global*3];
                        geoArray[sp].theta=angles[proj_global*3+1];
                        geoArray[sp].psi  =angles[proj_global*3+2];
                        Point3D source, deltaU, deltaV, uvOrigin;
                        computeDeltas_Siddon_parallel(geoArray[sp],geoArray[sp].alpha,proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                        params.data[4*j]=uvOrigin;
                        params.data[4*j+1]=deltaU;
                        params.data[4*j+2]=deltaV;
                        params.data[4*j+3]=source;
                    }
                    sycl::queue& q = tigreSyclQueue(gpuids[dev]);
                    float* detector = dProjection[dev];
                    const float* tex = dImage[dev];
                    const ProjParamsParallel p = params;
                    const int nangles_dev = nangles_device;
                    q.parallel_for(sycl::range<3>(PROJ_PER_BLOCK, geo.nDetecV, geo.nDetecU), [=](sycl::id<3> id){
                        kernelPixelDetector_parallel(geoSplit, detector, i, nangles_dev, tex, p.data, id);
                    });
                    float* dst = (sp==0) ? result[proj_start] : &partial[dev*num_pixels_block];
                    copies[dev] = q.memcpy(dst, detector, (size_t)projection_this_block[dev]*geo.nDetecV*geo.nDetecU*sizeof(float));
                }
                for (dev=0;dev<deviceCount;dev++){
                    if (!projection_this_block[dev])
                        continue;
                    copies[dev].wait_and_throw();
                    if (sp>0){
                        // Add the contribution of this slab of the image
                        unsigned int proj_start=i*PROJ_PER_BLOCK+dev*nangles_device;
                        float* acc = result[proj_start];
                        const float* add = &partial[dev*num_pixels_block];
                        size_t n = (size_t)projection_this_block[dev]*geo.nDetecV*geo.nDetecU;
                        for (size_t e=0; e<n; e++)
                            acc[e]+=add[e];
                    }
                }
            }
        }
        onError.dismiss();  // from here on the function releases everything itself
        for (dev = 0; dev < deviceCount; dev++){
            sycl::queue& q = tigreSyclQueue(gpuids[dev]);
            q.wait_and_throw();
            sycl::free(dImage[dev], q);
            sycl::free(dProjection[dev], q);
        }
        freeGeoArray(splits,geoArray);
    } catch (const sycl::exception& e) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection_parallel:SYCL", e.what());
    }
    return 0;
}

/* This code precomputes The location of the source and the Delta U and delta V (in the warped space)
 * to compute the locations of the x-rays. While it seems verbose and overly-optimized,
 * it does saves about 30% of each of the kernel calls. Thats something!
 **/
void computeDeltas_Siddon_parallel(Geometry geo, float angles,int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){
    Point3D S;
    
    S.x  =geo.DSO[i];   S.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       S.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD[i]-geo.DSO[i]);   P.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       P.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pu0.x=-(geo.DSD[i]-geo.DSO[i]);   Pu0.y= geo.dDetecU*(1-((float)geo.nDetecU/2)+0.5);       Pu0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pv0.x=-(geo.DSD[i]-geo.DSO[i]);   Pv0.y= geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       Pv0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-1);
    // Geometric trasnformations:
    P.x=0;Pu0.x=0;Pv0.x=0;
    
    // Roll pitch yaw
    rollPitchYaw(geo,i,&P);
    rollPitchYaw(geo,i,&Pu0);
    rollPitchYaw(geo,i,&Pv0);
    //Now lets translate the points where they should be:
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);
    Pu0.x=Pu0.x-(geo.DSD[i]-geo.DSO[i]);
    Pv0.x=Pv0.x-(geo.DSD[i]-geo.DSO[i]);

    S.x=0;
    // Roll pitch yaw
    rollPitchYaw(geo,i,&S);
    //Now lets translate the points where they should be:
    S.x=S.x+geo.DSO[i];

    //1: Offset detector
    
    //P.x
    P.y  =P.y  +geo.offDetecU[i];    P.z  =P.z  +geo.offDetecV[i];
    Pu0.y=Pu0.y+geo.offDetecU[i];    Pu0.z=Pu0.z+geo.offDetecV[i];
    Pv0.y=Pv0.y+geo.offDetecU[i];    Pv0.z=Pv0.z+geo.offDetecV[i];
    //S doesnt need to chagne
    
    
    //3: Rotate using full Euler angles (ZYZ convention: alpha, theta, psi)
    Point3D Pfinal = P;
    Point3D Pfinalu0 = Pu0;
    Point3D Pfinalv0 = Pv0;
    Point3D S2 = S;
    
    eulerZYZ(geo, &Pfinal);
    eulerZYZ(geo, &Pfinalu0);
    eulerZYZ(geo, &Pfinalv0);
    eulerZYZ(geo, &S2);
    
    //2: Offset image (instead of offseting image, -offset everything else)
    
    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    Pfinalu0.x=Pfinalu0.x-geo.offOrigX[i];   Pfinalu0.y=Pfinalu0.y-geo.offOrigY[i];   Pfinalu0.z=Pfinalu0.z-geo.offOrigZ[i];
    Pfinalv0.x=Pfinalv0.x-geo.offOrigX[i];   Pfinalv0.y=Pfinalv0.y-geo.offOrigY[i];   Pfinalv0.z=Pfinalv0.z-geo.offOrigZ[i];
    S2.x=S2.x-geo.offOrigX[i];               S2.y=S2.y-geo.offOrigY[i];               S2.z=S2.z-geo.offOrigZ[i];
    
    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2;
    Pfinalu0.x=Pfinalu0.x+geo.sVoxelX/2;    Pfinalu0.y=Pfinalu0.y+geo.sVoxelY/2;        Pfinalu0.z=Pfinalu0.z+geo.sVoxelZ/2;
    Pfinalv0.x=Pfinalv0.x+geo.sVoxelX/2;    Pfinalv0.y=Pfinalv0.y+geo.sVoxelY/2;        Pfinalv0.z=Pfinalv0.z+geo.sVoxelZ/2;
    S2.x      =S2.x+geo.sVoxelX/2;          S2.y      =S2.y+geo.sVoxelY/2;              S2.z      =S2.z      +geo.sVoxelZ/2;
    
    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S2.x      =S2.x/geo.dVoxelX;          S2.y      =S2.y/geo.dVoxelY;            S2.z      =S2.z/geo.dVoxelZ;
    
    
    
    //5. apply COR. Wherever everything was, now its offesetd by a bit
    float CORx, CORy;
    CORx=-geo.COR[i]*sin(geo.alpha)/geo.dVoxelX;
    CORy= geo.COR[i]*cos(geo.alpha)/geo.dVoxelY;
    Pfinal.x+=CORx;   Pfinal.y+=CORy;
    Pfinalu0.x+=CORx;   Pfinalu0.y+=CORy;
    Pfinalv0.x+=CORx;   Pfinalv0.y+=CORy;
    S2.x+=CORx; S2.y+=CORy;
    
    // return
    
    *uvorigin=Pfinal;
    
    deltaU->x=Pfinalu0.x-Pfinal.x;
    deltaU->y=Pfinalu0.y-Pfinal.y;
    deltaU->z=Pfinalu0.z-Pfinal.z;
    
    deltaV->x=Pfinalv0.x-Pfinal.x;
    deltaV->y=Pfinalv0.y-Pfinal.y;
    deltaV->z=Pfinalv0.z-Pfinal.z;
    
    *source=S2;
}

#ifndef PROJECTION_HPP

float maxDistanceCubeXY(Geometry geo, float alpha,int i){
    ///////////
    // Compute initial "t" so we access safely as less as out of bounds as possible.
    //////////
    
    
    float maxCubX,maxCubY;
    // Forgetting Z, compute max distance: diagonal+offset
    maxCubX=(geo.sVoxelX/2+ abs(geo.offOrigX[i]))/geo.dVoxelX;
    maxCubY=(geo.sVoxelY/2+ abs(geo.offOrigY[i]))/geo.dVoxelY;
    
    return geo.DSO[i]/geo.dVoxelX-sqrt(maxCubX*maxCubX+maxCubY*maxCubY);
    
}

#endif
