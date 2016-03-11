

#ifndef TYPES_CBCT
#define TYPES_CBCT
struct  Geometry {
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
};

 struct Point3D{
    float x;
    float y;
    float z;
};
#endif


#ifndef PROJECTION_HPP
#define PROJECTION_HPP

int interpolation_projection(float const * const img, Geometry geo, float** result,float const * const alphas,int nalpha);
float computeMaxLength(Geometry geo, float alpha);
void computeDeltas(Geometry geo, float alpha,int i, Point3D* uvorigin, Point3D* deltaU, Point3D* deltaV, Point3D* source);

float maxDistanceCubeXY(Geometry geo, float alpha,int i);

// below, not used
Geometry nomralizeGeometryImage(Geometry geo);
#endif