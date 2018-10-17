#ifndef TYPES
#define TYPES

typedef struct Geometry{
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
    float* dRoll;
    float* dPitch;
    float* dYaw;
    // The base unit we are working with in mm.
    float unitX;
    float unitY;
    float unitZ;

    //projection angle
    float alpha;
    // Centre of Rotation correction.
    float* COR;
    //Maximum length of cube
    float maxLength;
     //User option
    float accuracy;
}CBCT_Geometry;

typedef struct Point3D{
    float x;
    float y;
    float z;
}Point3D;

#endif