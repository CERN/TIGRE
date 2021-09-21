/*-------------------------------------------------------------------------
 *
 * Header CUDA function for backrpojection using matched weigts for CBCT
 *
 *
 * CODE by  Ander Biguri
 *          Optimized and modified by RB
---------------------------------------------------------------------------
---------------------------------------------------------------------------
Copyright (c) 2015, University of Bath and CERN- European Organization for 
Nuclear Research
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
 ---------------------------------------------------------------------------

Contact: tigre.toolbox@gmail.com
Codes  : https://github.com/CERN/TIGRE
--------------------------------------------------------------------------- 
 */




#include "voxel_backprojection.hpp"
#include "types_TIGRE.hpp"
#include "GpuIds.hpp"


#ifndef BACKPROJECTION2_HPP
#define BACKPROJECTION2_HPP

int voxel_backprojection2(float  *  projections, Geometry geo, float* result,float const * const alphas,int nalpha, const GpuIds& gpuids);
void computeDeltasCube(Geometry geo, float alpha,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S);
void splitCTbackprojection(const GpuIds& gpuids,Geometry geo,int nalpha, unsigned int* split_image, unsigned int * split_projections);
void computeDeltasCube(Geometry geo, int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S);
void createGeoArray(unsigned int image_splits, Geometry geo,Geometry* geoArray, unsigned int nangles);
void freeGeoArray(unsigned int splits,Geometry* geoArray);
void checkFreeMemory(const GpuIds& gpuids, size_t *mem_GPU_global);
#endif