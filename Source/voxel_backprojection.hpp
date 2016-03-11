
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

#ifndef BACKPROJECTION_HPP
#define BACKPROJECTION_HPP

int voxel_backprojection(float const * const projections, Geometry geo, float* result,float const * const alphas,int nalpha);
void computeDeltasCube(Geometry geo, float alpha,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ);
#endif