/*-------------------------------------------------------------------------
 *
 * Header CUDA function for backrpojection using matched weigts for CBCT
 *
 *
 * CODE by  Ander Biguri
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
#ifndef TYPES_CBCT
typedef struct {
    // Geometry assumptions:
    //  -> Origin is at (0,0,0). Image center is there +offOrig
    //  -> at angle 0, source + image centre (without the offset) + detector centre (without offset) 
    //     are aligned in the Y_Z plane.
    //  -> detector is orthonormal to projection plane.
    
    //Parameters part of the image geometry
    int   nVoxelX, nVoxelY, nVoxelZ;
    float sVoxelX, sVoxelY, sVoxelZ;
    float dVoxelX, dVoxelY, dVoxelZ;
    float *offOrigX,*offOrigY,*offOrigZ;
    float DSO;
    // Parameters  of the Detector.
    int   nDetecU, nDetecV;
    float sDetecU, sDetecV;
    float dDetecU, dDetecV;
    float *offDetecU, *offDetecV;
    float DSD;
    
    // The base unit we are working with in mm. 
    float unitX;
    float unitY;
    float unitZ;
    
    //projection angle
    float alpha;
    
    // Centre of Rotation correction.
    float COR;
    //Maximum length of cube
    float maxLength;
     //User option
    float accuracy;
}Geometry;

typedef struct{
    float x;
    float y;
    float z;
}Point3D;
#define TYPES_CBCT
#endif

#ifndef BACKPROJECTION2_HPP
#define BACKPROJECTION2_HPP

int voxel_backprojection2(float const * const projections, Geometry geo, float* result,float const * const alphas,int nalpha);
void computeDeltasCube(Geometry geo, float alpha,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ);
#endif