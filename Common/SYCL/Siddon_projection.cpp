/*-------------------------------------------------------------------------
 *
 * SYCL functions for ray-voxel intersection based projection
 *
 * This file has the necessary fucntiosn to perform X-ray CBCT projection
 * operation given a geaometry, angles and image. It usesthe so-called
 * Jacobs algorithm to compute efficiently the length of the x-rays over
 * voxel space.
 *
 * CODE by       Ander Biguri
 *               Sepideh Hatamikia (arbitrary rotation)
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
#include "TIGRE_common.hpp"
#include "tigre_sycl.hpp"
#include "tigre_release_on_error_sycl.hpp"
#include "page_locking_sycl.hpp"
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
struct ProjParams { Point3D data[4*PROJ_PER_BLOCK]; };

static inline void kernelPixelDetector(const Geometry& geo,
        float* detector,
        const int currProjSetNumber,
        const int totalNoOfProjections,
        const float* tex,
        const Point3D* projParamsArrayDev,
        sycl::id<3> id){
    
    
    unsigned long long u = id[2];
    unsigned long long v = id[1];
    unsigned long long projNumber=id[0];
    
    
    if (u>= geo.nDetecU || v>= geo.nDetecV || projNumber>=PROJ_PER_BLOCK)
        return;
    
#if IS_FOR_MATLAB_TIGRE
    size_t idx =  (size_t)(u * (unsigned long long)geo.nDetecV + v)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#else
    size_t idx =  (size_t)(v * (unsigned long long)geo.nDetecU + u)+ projNumber*(unsigned long long)geo.nDetecV *(unsigned long long)geo.nDetecU ;
#endif
    unsigned long indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array (for a given GPU)

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
    ///////
    // Siddon's ray-voxel intersection, optimized as in doi=10.1.1.55.7516
    //////
    // Also called Jacobs algorithms
    Point3D ray;
    // vector of Xray
    ray.x=pixel1D.x-source.x;
    ray.y=pixel1D.y-source.y;
    ray.z=pixel1D.z-source.z;
    float eps=0.001;
    ray.x=(sycl::fabs(ray.x)<eps)? 0 : ray.x;
    ray.y=(sycl::fabs(ray.y)<eps)? 0 : ray.y; 
    ray.z=(sycl::fabs(ray.z)<eps)? 0 : ray.z; 
    // This variables are ommited because
    // bx,by,bz ={0,0,0}
    // dx,dy,dz ={1,1,1}
    // compute parameter values for x-ray parametric equation. eq(3-10)
    float axm,aym,azm;
    float axM,ayM,azM;
    // In the paper Nx= number of X planes-> Nvoxel+1
   
    axm=sycl::fmin(((-source.x)/(ray.x)),((geo.nVoxelX-source.x)/(ray.x)));
    aym=sycl::fmin(((-source.y)/(ray.y)),((geo.nVoxelY-source.y)/(ray.y)));
    azm=sycl::fmin(((-source.z)/(ray.z)),((geo.nVoxelZ-source.z)/(ray.z)));
    axM=sycl::fmax(((-source.x)/(ray.x)),((geo.nVoxelX-source.x)/(ray.x)));
    ayM=sycl::fmax(((-source.y)/(ray.y)),((geo.nVoxelY-source.y)/(ray.y)));
    azM=sycl::fmax(((-source.z)/(ray.z)),((geo.nVoxelZ-source.z)/(ray.z)));
    
    float am=sycl::fmax(sycl::fmax(axm,aym),azm);
    float aM=sycl::fmin(sycl::fmin(axM,ayM),azM);
    
    // line intersects voxel space ->   am<aM
    // The ray never enters the volume. RETURN HERE. Without it, execution
    // falls into the index block below, whose imin..kmax are only meaningful
    // for a ray that does enter. For a ray the eps guard has snapped to run
    // exactly along one voxel axis, with the source OUTSIDE the volume on
    // one of its zero axes, the two quotients for that axis share a sign, so
    // fminf and fmaxf collapse (axm == axM == -inf) and aM inherits -inf.
    // The other zero component then gives aM*ray = (-inf)*0 = NaN, Np is
    // NaN, and (unsigned long)NaN is 2^63 on CUDA: that thread iterates
    // ~9e18 times and the kernel never returns -- 100% SM, idle memory
    // controller, no error, no result. Measured: 2 rays in 261 million on an
    // odd x odd detector with an off-axis object; see
    // Python/tests/test_siddon_missing_ray.py and
    // Frontispiece/siddon_missed_ray_hang.png. A NEGATIVE Np does not do
    // this -- CUDA saturates negatives and -inf to 0, a zero-trip loop.
    if (am>=aM){
        detector[idx]=0;
        return;
    }

    // Walk the voxels the ray crosses, from where it enters the box (am) to where
    // it leaves (aM), adding up length x value (#859, as in the CUDA kernel).
    //
    // The first voxel is the one the ray is in just after entering. When the
    // entry point lies on a voxel boundary of another axis too, the direction
    // decides the side: moving up, floor(); moving down, ceil()-1. The next
    // crossing on each axis is then computed from that voxel's own faces, so the
    // two always agree. Before, the first voxel came from a midpoint and the
    // first crossings from ceil/floor of the entry point, separately; on such
    // rays they disagreed, giving a wrong first voxel or the whole ray one voxel
    // off along an axis. Rare for the whole image, but common when the image is
    // split in pieces, where many rays enter a piece through its z faces.
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

    float maxlength=sycl::sqrt(ray.x*ray.x*geo.dVoxelX*geo.dVoxelX+ray.y*ray.y*geo.dVoxelY*geo.dVoxelY+ray.z*ray.z*geo.dVoxelZ*geo.dVoxelZ);
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
    detector[idx]=sum*maxlength;
}


int siddon_ray_projection(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // GPU arrays (#857): the pointers would be CUDA device memory, which a SYCL device
    // cannot read. tigre.Ax and tigre.Atb pass such arrays through the host instead.
    if (gpuData) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection:GPUdata","GPU arrays (PyTorch CUDA tensors, CuPy arrays) are not supported by the SYCL build: pass host (NumPy) arrays.\n");
    }
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection:GPUselect","There are no available SYCL device(s)\n");
    }
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Ax:Siddon_projection:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
    }
    int dev;
    // If an error raises (Python build) before the end of this function, onError
    // releases what it holds; see tigre_release_on_error_sycl.hpp.
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
                tigreNotEnoughGpuMemory("Ax:Siddon_projection:GPUmemory","Forward projection (Siddon)",4*PROJ_PER_BLOCK*mem_proj+mem_slice,mem_GPU_global);
            size_t mem_free=mem_GPU_global-4*PROJ_PER_BLOCK*mem_proj;
            splits=mem_image/mem_free+1;// Ceil of the truncation
            if (splits>geo.nVoxelZ) splits=geo.nVoxelZ;  // at least one slice per piece
        }
        // Allocate the image (a piece of it if split) and one block of projections per device. If that fails
        // although the estimate above fits, split the image in more pieces and try again, down to one slice per
        // piece, as the CUDA code does.
        size_t num_pixels_block = (size_t)PROJ_PER_BLOCK*(size_t)geo.nDetecU*(size_t)geo.nDetecV;
        std::vector<int> devs(deviceCount);
        for (dev = 0; dev < deviceCount; dev++) devs[dev] = gpuids[dev];
        std::vector<float*> dImage(deviceCount), dProjection(deviceCount);
        Geometry* geoArray;
        while (true) {
            // As many pieces as their size needs: no empty pieces at the end
            splits=(geo.nVoxelZ+(geo.nVoxelZ+splits-1)/splits-1)/((geo.nVoxelZ+splits-1)/splits);
            geoArray = (Geometry*)malloc(splits*sizeof(Geometry));
            splitImage(splits,geo,geoArray,nangles);
            size_t num_voxels_split = (size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
            if (tigre::tryMallocPieces(devs.data(), deviceCount, num_voxels_split, num_pixels_block, dImage.data(), dProjection.data()))
                break;
            const unsigned int piece = geoArray[0].nVoxelZ;
            freeGeoArray(splits, geoArray);
            if (piece<=1 || splits>=geo.nVoxelZ)
                tigreNotEnoughGpuMemory("Ax:Siddon_projection:GPUmemory","Forward projection (Siddon)",(num_pixels_block+(size_t)geo.nVoxelX*geo.nVoxelY)*sizeof(float),tigreSyclFreeMemory(gpuids[0]));
            while (splits<geo.nVoxelZ && (geo.nVoxelZ+splits-1)/splits>=piece) splits++;
        }
        onError.add([geoArray, splits]{ freeGeoArray(splits, geoArray); });
        for (dev = 0; dev < deviceCount; dev++){
            sycl::queue* q = &tigreSyclQueue(gpuids[dev]);
            onError.addQueue(q);
            float* image_d = dImage[dev];
            float* projection_d = dProjection[dev];
            onError.add([image_d, projection_d, q]{ sycl::free(image_d, *q); sycl::free(projection_d, *q); });
        }

        // Page-lock the image for the copies to the devices only when it is worth it (split, or more than 2 devices)
        bool img_registered = false;
        if (tigre::worthPageLocking(tigre::hostRegisterSupported(devs.data(), deviceCount), gpuData, splits, deviceCount)) {
            img_registered = tigre::hostRegister(img, mem_image, devs.data(), deviceCount);
            if (img_registered) onError.add([img, devs]{ tigre::hostUnregister(img, devs.data(), (int)devs.size()); });
        }
        std::vector<float> partial(splits>1 ? num_pixels_block*deviceCount : 0);

        // Split the angles between devices
        int nangles_device=(nangles+deviceCount-1)/deviceCount;
        unsigned int noOfKernelCalls = (nangles_device+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_BLOCK

        ProjParams params;
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
                        computeDeltas_Siddon(geoArray[sp],proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                        params.data[4*j]=uvOrigin;
                        params.data[4*j+1]=deltaU;
                        params.data[4*j+2]=deltaV;
                        params.data[4*j+3]=source;
                    }
                    sycl::queue& q = tigreSyclQueue(gpuids[dev]);
                    float* detector = dProjection[dev];
                    const float* tex = dImage[dev];
                    const ProjParams p = params;
                    const int nangles_dev = nangles_device;
                    q.parallel_for(sycl::range<3>(PROJ_PER_BLOCK, geo.nDetecV, geo.nDetecU), [=](sycl::id<3> id){
                        kernelPixelDetector(geoSplit, detector, i, nangles_dev, tex, p.data, id);
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
        if (img_registered) tigre::hostUnregister(img, devs.data(), deviceCount);
    } catch (const sycl::exception& e) {
        mexErrMsgIdAndTxt("Ax:Siddon_projection:SYCL", e.what());
    }
    return 0;
}


/* This code generates the geometries needed to split the image properly in
 * cases where the entire image does not fit in the memory of the GPU
 **/
void splitImage(unsigned int splits,Geometry geo,Geometry* geoArray, unsigned int nangles){
    
    unsigned long splitsize=(geo.nVoxelZ+splits-1)/splits;// ceil if not divisible
    for(unsigned int sp=0;sp<splits;sp++){
        geoArray[sp]=geo;
        // All of them are splitsize, but the last one; signed, so that a piece beyond the image is empty
        const long long rest=(long long)geo.nVoxelZ-(long long)splitsize*sp;
        geoArray[sp].nVoxelZ=(int)(rest<=0 ? 0 : (rest<(long long)splitsize ? rest : (long long)splitsize));
        geoArray[sp].sVoxelZ= geoArray[sp].nVoxelZ* geoArray[sp].dVoxelZ;
        
        // We need to redefine the offsets, as now each subimage is not aligned in the origin.
        geoArray[sp].offOrigZ=(float *)malloc(nangles*sizeof(float));
        for (unsigned int i=0;i<nangles;i++){
            geoArray[sp].offOrigZ[i]=geo.offOrigZ[i]-geo.sVoxelZ/2+sp*geoArray[0].sVoxelZ+geoArray[sp].sVoxelZ/2;
        }
        
    }
    
}

/* This code precomputes The location of the source and the Delta U and delta V (in the warped space)
 * to compute the locations of the x-rays. While it seems verbose and overly-optimized,
 * it does saves about 30% of each of the kernel calls. Thats something!
 **/
void computeDeltas_Siddon(Geometry geo,int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){

    
    Point3D S;
    S.x=geo.DSO[i];
    S.y=0;
    S.z=0;
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD[i]-geo.DSO[i]);   P.y  = geo.dDetecU*(-((double)geo.nDetecU/2.0)+0.5);       P.z  = geo.dDetecV*(((double)geo.nDetecV/2.0)-0.5);
    Pu0.x=0;                          Pu0.y= geo.dDetecU;                                    Pu0.z= 0;
    Pv0.x=0;                          Pv0.y= 0;                                              Pv0.z= geo.dDetecV*(-1);

    // Geometric transformations:
    // Now we have the Real world (OXYZ) coordinates of the bottom corner and its two neighbours.
    // The objective is to get a position of the detector in a coordinate system where:
    // 1-units are voxel size (in each direction can be different)
    // 2-The image has the its first voxel at (0,0,0)
    // 3-The image never rotates

    // To do that, we need to compute the "deltas" the detector, or "by how much
    // (in new xyz) does the voxels change when and index is added". To do that
    // several geometric steps needs to be changed

    //1.Roll,pitch,jaw
    // The detector can have a small rotation.
    // according to
    //"A geometric calibration method for cone beam CT systems" Yang K1, Kwan AL, Miller DF, Boone JM. Med Phys. 2006 Jun;33(6):1695-706.
    // Only the Z rotation will have a big influence in the image quality when they are small.
    // Still all rotations are supported

    // To roll pitch jaw, the detector has to be in centered in OXYZ.
    // NB: do not apply offsets to Pu0 and Pv0: they are directions, and are invariant through translations
    P.x=0;

    // Roll pitch yaw
    rollPitchYaw(geo,i,&P);
    rollPitchYaw(geo,i,&Pu0);
    rollPitchYaw(geo,i,&Pv0);
    //Now let's translate the points where they should be:
    // NB: do not apply offsets to Pu0 and Pv0: they are directions, and are invariant through translations
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);

    //1: Offset detector


    //S doesnt need to chagne


    //3: Rotate (around z)!
    Point3D Pfinal, Pfinalu0, Pfinalv0;
    Pfinal.x  =P.x;
    Pfinal.y  =P.y  +geo.offDetecU[i]; Pfinal.z  =P.z  +geo.offDetecV[i];
    Pfinalu0 = Pu0;
    Pfinalv0 = Pv0;

    eulerZYZ(geo,&Pfinal);
    eulerZYZ(geo,&Pfinalu0);
    eulerZYZ(geo,&Pfinalv0);
    eulerZYZ(geo,&S);

    //2: Offset image (instead of offseting image, -offset everything else)
    // NB: do not apply offsets to Pfinalu0 and Pfinalv0: they are directions, and are invariant through translations

    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    S.x=S.x-geo.offOrigX[i];               S.y=S.y-geo.offOrigY[i];               S.z=S.z-geo.offOrigZ[i];

    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2;
    S.x      =S.x+geo.sVoxelX/2;          S.y      =S.y+geo.sVoxelY/2;              S.z      =S.z      +geo.sVoxelZ/2;

    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S.x      =S.x/geo.dVoxelX;          S.y      =S.y/geo.dVoxelY;            S.z      =S.z/geo.dVoxelZ;


    //mexPrintf("COR: %f \n",geo.COR[i]);
    //5. apply COR. Wherever everything was, now its offesetd by a bit
    // NB: do not apply offsets to Pfinalu0 and Pfinalv0: they are directions, and are invariant through translations
    double CORx, CORy;
    CORx=-geo.COR[i]*sin(geo.alpha)/geo.dVoxelX;
    CORy= geo.COR[i]*cos(geo.alpha)/geo.dVoxelY;
    Pfinal.x+=CORx;   Pfinal.y+=CORy;
    S.x+=CORx; S.y+=CORy;

    // return

    *uvorigin=Pfinal;

    *deltaU=Pfinalu0;
    *deltaV=Pfinalv0;
    
    *source=S;
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
void rollPitchYaw(Geometry geo,int i, Point3D* point){
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x=cos(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.x
            +(cos(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) - sin(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.y
            +(cos(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) + sin(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.z;
    
    point->y=sin(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.x
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) + cos(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.y
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) - cos(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.z;
    
    point->z=-sin(geo.dPitch[i])*auxPoint.x
            +cos(geo.dPitch[i])*sin(geo.dYaw[i])*auxPoint.y
            +cos(geo.dPitch[i])*cos(geo.dYaw[i])*auxPoint.z;
    
}
void eulerZYZ(Geometry geo, Point3D* point){
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x=(+cos(geo.alpha)*cos(geo.theta)*cos(geo.psi)-sin(geo.alpha)*sin(geo.psi))*auxPoint.x+
            (-cos(geo.alpha)*cos(geo.theta)*sin(geo.psi)-sin(geo.alpha)*cos(geo.psi))*auxPoint.y+
            cos(geo.alpha)*sin(geo.theta)*auxPoint.z;
    
    point->y=(+sin(geo.alpha)*cos(geo.theta)*cos(geo.psi)+cos(geo.alpha)*sin(geo.psi))*auxPoint.x+
            (-sin(geo.alpha)*cos(geo.theta)*sin(geo.psi)+cos(geo.alpha)*cos(geo.psi))*auxPoint.y+
            sin(geo.alpha)*sin(geo.theta)*auxPoint.z;
    
    point->z=-sin(geo.theta)*cos(geo.psi)*auxPoint.x+
            sin(geo.theta)*sin(geo.psi)*auxPoint.y+
            cos(geo.theta)*auxPoint.z;
    
    
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
void checkFreeMemory(const GpuIds& gpuids, size_t *mem_GPU_global){
    // Free memory where the SYCL runtime reports it, otherwise the global memory
    // of the device (95%); never more than the largest single allocation.
    const int deviceCount = gpuids.GetLength();
    for (int dev = 0; dev < deviceCount; dev++){
        size_t mem = tigreSyclCheckFreeMemory(gpuids[dev],"Ax:Siddon_projection:GPUmemory");
        if(dev==0) *mem_GPU_global=mem;
        *mem_GPU_global=(mem<*mem_GPU_global)?mem:*mem_GPU_global;
    }
}
#endif
