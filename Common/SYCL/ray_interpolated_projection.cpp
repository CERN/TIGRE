/*-------------------------------------------------------------------------
 *
 * SYCL functions for texture-memory interpolation based projection
 *
 * This file has the necessary fucntiosn to perform X-ray CBCT projection
 * operation given a geaometry, angles and image. It uses the 3D texture
 * memory linear interpolation to uniformily sample a path to integrate the
 * X-rays.
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
#include "ray_interpolated_projection.hpp"
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
struct ProjParamsInterp { Point3D data[4*PROJ_PER_BLOCK]; float floats[2*PROJ_PER_BLOCK]; };

template<bool sphericalrotation>
static inline void kernelPixelDetector(const Geometry& geo,
        float* detector,
        const int currProjSetNumber,
        const int totalNoOfProjections,
        const float* tex,
        const bool linearInterpolation,
        const Point3D* projParamsArrayDev,
        const float* projFloatsArrayDev,
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

    unsigned long indAlpha = currProjSetNumber*PROJ_PER_BLOCK+projNumber;  // This is the ABSOLUTE projection number in the projection array
    
    if(indAlpha>=totalNoOfProjections)
        return;
    
    Point3D uvOrigin = projParamsArrayDev[4*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
    Point3D deltaU = projParamsArrayDev[4*projNumber+1];
    Point3D deltaV = projParamsArrayDev[4*projNumber+2];
    Point3D source = projParamsArrayDev[4*projNumber+3];
    
    float DSO = projFloatsArrayDev[2*projNumber+0];
    float cropdist_init = projFloatsArrayDev[2*projNumber+1];
    
    
    
    /////// Get coordinates XYZ of pixel UV
    unsigned long pixelV = geo.nDetecV-v-1;
    unsigned long pixelU = u;
    
    
    float vectX,vectY,vectZ;
    Point3D P;
    P.x=(uvOrigin.x+pixelU*deltaU.x+pixelV*deltaV.x);
    P.y=(uvOrigin.y+pixelU*deltaU.y+pixelV*deltaV.y);
    P.z=(uvOrigin.z+pixelU*deltaU.z+pixelV*deltaV.z);
    
    // Length is the ray length in normalized space
    float length=sycl::sqrt((source.x-P.x)*(source.x-P.x)+(source.y-P.y)*(source.y-P.y)+(source.z-P.z)*(source.z-P.z));
    //now legth is an integer of Nsamples that are required on this line
    length=sycl::ceil(((length)/(geo.accuracy)));//Divide the directional vector by an integer
    vectX=((P.x -source.x)/(length));
    vectY=((P.y -source.y)/(length));
    vectZ=((P.z -source.z)/(length));
    
    
//     //Integrate over the line
    float tx,ty,tz;
    float sum=0;
    float i;
    
    
    
//  Because I have no idea how to efficiently cutoff the legth path in 3D, a very upper limit is computed (see maxdistanceCuboid)
//  for the 3D case. However it would be bad to lose performance in the 3D case
//  TODO: can ge really improve this?
    if (sphericalrotation){
        if ((2*DSO/sycl::fmin(sycl::fmin(geo.dVoxelX,geo.dVoxelY),geo.dVoxelZ)+cropdist_init)/geo.accuracy  <   length)
            length=sycl::ceil((2*DSO/sycl::fmin(sycl::fmin(geo.dVoxelX,geo.dVoxelY),geo.dVoxelZ)+cropdist_init)/geo.accuracy);
    }
    else{
        if ((2*DSO/sycl::fmin(geo.dVoxelX,geo.dVoxelY)+cropdist_init)/geo.accuracy  <   length)
            length=sycl::ceil((2*DSO/sycl::fmin(geo.dVoxelX,geo.dVoxelY)+cropdist_init)/geo.accuracy);
    }
    
    
    //Length is not actually a length, but the amount of memreads with given accuracy ("samples per voxel")
    for (i=sycl::floor(cropdist_init/geo.accuracy); i<=length; i=i+1){
        tx=vectX*i+source.x;
        ty=vectY*i+source.y;
        tz=vectZ*i+source.z;
        
        sum += linearInterpolation ? tigreFetchLinear3D(tex, geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ, tx+0.5f, ty+0.5f, tz+0.5f)
                                   : tigreFetchPoint3D (tex, geo.nVoxelX, geo.nVoxelY, geo.nVoxelZ, tx+0.5f, ty+0.5f, tz+0.5f); // this line is 94% of time.
    }
    
    float deltalength=sycl::sqrt((vectX*geo.dVoxelX)*(vectX*geo.dVoxelX)+
            (vectY*geo.dVoxelY)*(vectY*geo.dVoxelY)+
            (vectZ*geo.dVoxelZ)*(vectZ*geo.dVoxelZ) );
    
    detector[idx]=sum*deltalength;
}


int interpolation_projection(float* img, Geometry geo, float** result,float const * const angles,int nangles, const GpuIds& gpuids, bool gpuData){
    // GPU arrays (#857): the pointers would be CUDA device memory, which a SYCL device
    // cannot read. tigre.Ax and tigre.Atb pass such arrays through the host instead.
    if (gpuData) {
        mexErrMsgIdAndTxt("Ax:Interpolated_projection:GPUdata","GPU arrays (PyTorch CUDA tensors, CuPy arrays) are not supported by the SYCL build: pass host (NumPy) arrays.\n");
    }
    // Prepare for MultiGPU
    int deviceCount = gpuids.GetLength();
    if (deviceCount == 0) {
        mexErrMsgIdAndTxt("Ax:Interpolated_projection:GPUselect","There are no available SYCL device(s)\n");
    }
    if (!gpuids.AreEqualDevices()) {
        mexWarnMsgIdAndTxt("Ax:Interpolated_projection:GPUselect","Detected one (or more) different GPUs.\n This code is not smart enough to separate the memory GPU wise if they have different computational times or memory limits.\n First GPU parameters used. If the code errors you might need to change the way GPU selection is performed.");
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
                tigreNotEnoughGpuMemory("Ax:Interpolated_projection:GPUmemory","Forward projection (interpolated)",4*PROJ_PER_BLOCK*mem_proj+mem_slice,mem_GPU_global);
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
        splitImageInterp(splits,geo,geoArray,nangles);
        onError.add([geoArray, splits]{ freeGeoArray(splits, geoArray); });

        // Allocate the image (a slab of it if split) and one block of projections per device
        size_t num_voxels_split = (size_t)geoArray[0].nVoxelX*(size_t)geoArray[0].nVoxelY*(size_t)geoArray[0].nVoxelZ;
        size_t num_pixels_block = (size_t)PROJ_PER_BLOCK*(size_t)geo.nDetecU*(size_t)geo.nDetecV;
        std::vector<float*> dImage(deviceCount), dProjection(deviceCount);
        for (dev = 0; dev < deviceCount; dev++){
            sycl::queue* q = &tigreSyclQueue(gpuids[dev]);
            onError.addQueue(q);
            float* image_d = tigre::deviceMalloc<float>(num_voxels_split, gpuids[dev], "Ax:Interpolated_projection:malloc", "Forward projection (interpolated)");
            onError.add([image_d, q]{ sycl::free(image_d, *q); });
            float* projection_d = tigre::deviceMalloc<float>(num_pixels_block, gpuids[dev], "Ax:Interpolated_projection:malloc", "Forward projection (interpolated)");
            onError.add([projection_d, q]{ sycl::free(projection_d, *q); });
            dImage[dev] = image_d;
            dProjection[dev] = projection_d;
        }
        std::vector<float> partial(splits>1 ? num_pixels_block*deviceCount : 0);

        // Split the angles between devices
        int nangles_device=(nangles+deviceCount-1)/deviceCount;
        unsigned int noOfKernelCalls = (nangles_device+PROJ_PER_BLOCK-1)/PROJ_PER_BLOCK;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_BLOCK

        ProjParamsInterp params;
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
                    float is_spherical=0;
                    for(unsigned int j=0; j<PROJ_PER_BLOCK; j++){
                        unsigned int proj_global=(i*PROJ_PER_BLOCK+j)+dev*nangles_device;
                        if (proj_global>=nangles)
                            break;
                        if ((i*PROJ_PER_BLOCK+j)>=nangles_device)
                            break;
                        geoArray[sp].alpha=angles[proj_global*3];
                        geoArray[sp].theta=angles[proj_global*3+1];
                        geoArray[sp].psi  =angles[proj_global*3+2];
                        is_spherical+=std::abs(geoArray[sp].theta)+std::abs(geoArray[sp].psi);
                        //precomute distances for faster execution
                        float maxdist=maxdistanceCuboid(geoArray[sp],proj_global);
                        //Precompute per angle constant stuff for speed
                        Point3D source, deltaU, deltaV, uvOrigin;
                        computeDeltas(geoArray[sp], proj_global, &uvOrigin, &deltaU, &deltaV, &source);
                        //Ray tracing!
                        params.data[4*j]=uvOrigin;   // 6*j because we have 6 Point3D values per projection
                        params.data[4*j+1]=deltaU;
                        params.data[4*j+2]=deltaV;
                        params.data[4*j+3]=source;
                        params.floats[2*j]=geo.DSO[proj_global];
                        params.floats[2*j+1]=floor(maxdist);
                    }
                    sycl::queue& q = tigreSyclQueue(gpuids[dev]);
                    float* detector = dProjection[dev];
                    const float* tex = dImage[dev];
                    const ProjParamsInterp p = params;
                    const int nangles_dev = nangles_device;
                    // accuracy>1 uses point sampling, as the CUDA texture setup does
                    const bool linear = !(geo.accuracy>1);
                    if (!is_spherical){
                        q.parallel_for(sycl::range<3>(PROJ_PER_BLOCK, geo.nDetecV, geo.nDetecU), [=](sycl::id<3> id){
                            kernelPixelDetector<false>(geoSplit, detector, i, nangles_dev, tex, linear, p.data, p.floats, id);
                        });
                    }
                    else{
                        q.parallel_for(sycl::range<3>(PROJ_PER_BLOCK, geo.nDetecV, geo.nDetecU), [=](sycl::id<3> id){
                            kernelPixelDetector<true>(geoSplit, detector, i, nangles_dev, tex, linear, p.data, p.floats, id);
                        });
                    }
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
        mexErrMsgIdAndTxt("Ax:Interpolated_projection:SYCL", e.what());
    }
    return 0;
}

/* This code generates the geometries needed to split the image properly in
 * cases where the entire image does not fit in the memory of the GPU
 **/
void splitImageInterp(unsigned int splits,Geometry geo,Geometry* geoArray, unsigned int nangles){
    
    unsigned long splitsize=(geo.nVoxelZ+splits-1)/splits;// ceil if not divisible
    for(unsigned int sp=0;sp<splits;sp++){
        geoArray[sp]=geo;
        // All of them are splitsize, but the last one, possible
        geoArray[sp].nVoxelZ=((sp+1)*splitsize<geo.nVoxelZ)?  splitsize:  geo.nVoxelZ-splitsize*sp;
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
void computeDeltas(Geometry geo,unsigned int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source){
    Point3D S;
    S.x=geo.DSO[i];
    S.y=0;
    S.z=0;
    
    //End point
    Point3D P,Pu0,Pv0;
    
    P.x  =-(geo.DSD[i]-geo.DSO[i]);   P.y  = geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       P.z  = geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pu0.x=-(geo.DSD[i]-geo.DSO[i]);   Pu0.y= geo.dDetecU*(1-((float)geo.nDetecU/2)+0.5);       Pu0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-0);
    Pv0.x=-(geo.DSD[i]-geo.DSO[i]);   Pv0.y= geo.dDetecU*(0-((float)geo.nDetecU/2)+0.5);       Pv0.z= geo.dDetecV*(((float)geo.nDetecV/2)-0.5-1);
    // Geomtric trasnformations:
    
    
    // Now we have the Real world (OXYZ) coordinates of the bottom corner and its two neighbours.
    // The obkjective is to get a position of the detector in a coordinate system where:
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
    P.x=0;Pu0.x=0;Pv0.x=0;
    
    // Roll pitch yaw
    rollPitchYaw(geo,i,&P);
    rollPitchYaw(geo,i,&Pu0);
    rollPitchYaw(geo,i,&Pv0);
    //Now ltes translate the detector coordinates to DOD (original position on real coordinate system:
    P.x=P.x-(geo.DSD[i]-geo.DSO[i]);
    Pu0.x=Pu0.x-(geo.DSD[i]-geo.DSO[i]);
    Pv0.x=Pv0.x-(geo.DSD[i]-geo.DSO[i]);
    //2: Offset detector
    
    
    //S doesnt need to chagne
    
    
    
    //3: Rotate around RZ RY RZ
    Point3D Pfinal, Pfinalu0, Pfinalv0;
    Pfinal.x  =P.x;
    Pfinal.y  =P.y  +geo.offDetecU[i]; Pfinal.z  =P.z  +geo.offDetecV[i];
    Pfinalu0.x=Pu0.x;
    Pfinalu0.y=Pu0.y  +geo.offDetecU[i]; Pfinalu0.z  =Pu0.z  +geo.offDetecV[i];
    Pfinalv0.x=Pv0.x;
    Pfinalv0.y=Pv0.y  +geo.offDetecU[i]; Pfinalv0.z  =Pv0.z  +geo.offDetecV[i];
    
    eulerZYZ(geo,&Pfinal);
    eulerZYZ(geo,&Pfinalu0);
    eulerZYZ(geo,&Pfinalv0);
    eulerZYZ(geo,&S);
    
    
    //3: Offset image (instead of offseting image, -offset everything else)
    
    Pfinal.x  =Pfinal.x-geo.offOrigX[i];     Pfinal.y  =Pfinal.y-geo.offOrigY[i];     Pfinal.z  =Pfinal.z-geo.offOrigZ[i];
    Pfinalu0.x=Pfinalu0.x-geo.offOrigX[i];   Pfinalu0.y=Pfinalu0.y-geo.offOrigY[i];   Pfinalu0.z=Pfinalu0.z-geo.offOrigZ[i];
    Pfinalv0.x=Pfinalv0.x-geo.offOrigX[i];   Pfinalv0.y=Pfinalv0.y-geo.offOrigY[i];   Pfinalv0.z=Pfinalv0.z-geo.offOrigZ[i];
    S.x=S.x-geo.offOrigX[i];                 S.y=S.y-geo.offOrigY[i];                 S.z=S.z-geo.offOrigZ[i];
    
    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    Pfinal.x  =Pfinal.x+geo.sVoxelX/2-geo.dVoxelX/2;      Pfinal.y  =Pfinal.y+geo.sVoxelY/2-geo.dVoxelY/2;          Pfinal.z  =Pfinal.z  +geo.sVoxelZ/2-geo.dVoxelZ/2;
    Pfinalu0.x=Pfinalu0.x+geo.sVoxelX/2-geo.dVoxelX/2;    Pfinalu0.y=Pfinalu0.y+geo.sVoxelY/2-geo.dVoxelY/2;        Pfinalu0.z=Pfinalu0.z+geo.sVoxelZ/2-geo.dVoxelZ/2;
    Pfinalv0.x=Pfinalv0.x+geo.sVoxelX/2-geo.dVoxelX/2;    Pfinalv0.y=Pfinalv0.y+geo.sVoxelY/2-geo.dVoxelY/2;        Pfinalv0.z=Pfinalv0.z+geo.sVoxelZ/2-geo.dVoxelZ/2;
    S.x       =S.x+geo.sVoxelX/2-geo.dVoxelX/2;           S.y       =S.y+geo.sVoxelY/2-geo.dVoxelY/2;               S.z       =S.z      +geo.sVoxelZ/2-geo.dVoxelZ/2;
    
    //4. Scale everything so dVoxel==1
    Pfinal.x  =Pfinal.x/geo.dVoxelX;      Pfinal.y  =Pfinal.y/geo.dVoxelY;        Pfinal.z  =Pfinal.z/geo.dVoxelZ;
    Pfinalu0.x=Pfinalu0.x/geo.dVoxelX;    Pfinalu0.y=Pfinalu0.y/geo.dVoxelY;      Pfinalu0.z=Pfinalu0.z/geo.dVoxelZ;
    Pfinalv0.x=Pfinalv0.x/geo.dVoxelX;    Pfinalv0.y=Pfinalv0.y/geo.dVoxelY;      Pfinalv0.z=Pfinalv0.z/geo.dVoxelZ;
    S.x       =S.x/geo.dVoxelX;           S.y       =S.y/geo.dVoxelY;             S.z       =S.z/geo.dVoxelZ;
    
    
    //mexPrintf("COR: %f \n",geo.COR[i]);
    //5. apply COR. Wherever everything was, now its offesetd by a bit.
//     Only wors for standard rotaiton, not aribtary axis rotation.
    float CORx, CORy;
    CORx=-geo.COR[i]*sin(geo.alpha)/geo.dVoxelX;
    CORy= geo.COR[i]*cos(geo.alpha)/geo.dVoxelY;
    Pfinal.x+=CORx;   Pfinal.y+=CORy;
    Pfinalu0.x+=CORx;   Pfinalu0.y+=CORy;
    Pfinalv0.x+=CORx;   Pfinalv0.y+=CORy;
    S.x+=CORx; S.y+=CORy;
    
    // return
    
    *uvorigin=Pfinal;
    
    deltaU->x=Pfinalu0.x-Pfinal.x;
    deltaU->y=Pfinalu0.y-Pfinal.y;
    deltaU->z=Pfinalu0.z-Pfinal.z;
    
    deltaV->x=Pfinalv0.x-Pfinal.x;
    deltaV->y=Pfinalv0.y-Pfinal.y;
    deltaV->z=Pfinalv0.z-Pfinal.z;
    
    *source=S;
}

float maxdistanceCuboid(Geometry geo,unsigned int i){
    ///////////
    // Compute initial "t" so we access safely as less as out of bounds as possible.
    //////////
    
    
    float maxCubX,maxCubY,maxCubZ;
    // Forgetting Z, compute mas distance: diagonal+offset
    maxCubX=(geo.nVoxelX/2+ abs(geo.offOrigX[i])/geo.dVoxelX);
    maxCubY=(geo.nVoxelY/2+ abs(geo.offOrigY[i])/geo.dVoxelY);
    maxCubZ=(geo.nVoxelZ/2+ abs(geo.offOrigZ[i])/geo.dVoxelZ);
    
    float a,b;
    a=geo.DSO[i]/geo.dVoxelX;
    b=geo.DSO[i]/geo.dVoxelY;
    
//  As the return of this value is in "voxel space", the source may have an elliptical curve.
//  The distance returned is the safe distance that can be skipped for a given angle alpha, before we need to start sampling.
    
    if (geo.theta==0.0f & geo.psi==0.0f) // Special case, it will make the code faster
        return std::max(a*b/sqrt(a*a*sin(geo.alpha)*sin(geo.alpha)+b*b*cos(geo.alpha)*cos(geo.alpha))-
                sqrt(maxCubX*maxCubX+maxCubY*maxCubY),0.0f);
    //TODO: think of more special cases?
    return std::max(geo.DSO[i]/std::max(std::max(geo.dVoxelX,geo.dVoxelY),geo.dVoxelZ)-sqrt(maxCubX*maxCubX+maxCubY*maxCubY+maxCubZ*maxCubZ),0.0f);

}
void rollPitchYaw(Geometry geo,unsigned int i, Point3D* point){
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
void eulerZYZ(Geometry geo,  Point3D* point){
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
//______________________________________________________________________________
//
//      Function:       checkFreeMemory
//
//      Description:    check available memory on devices
//______________________________________________________________________________
void checkFreeMemory(const GpuIds& gpuids, size_t *mem_GPU_global){
    // Free memory where the SYCL runtime reports it, otherwise the global memory
    // of the device (95%); never more than the largest single allocation.
    const int deviceCount = gpuids.GetLength();
    for (int dev = 0; dev < deviceCount; dev++){
        size_t mem = tigreSyclCheckFreeMemory(gpuids[dev],"Ax:Interpolated_projection:GPUmemory");
        if(dev==0) *mem_GPU_global=mem;
        *mem_GPU_global=(mem<*mem_GPU_global)?mem:*mem_GPU_global;
    }
}
