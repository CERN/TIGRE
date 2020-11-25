from libc.stdlib cimport malloc, free
cdef extern from "types_TIGRE.hpp":
    ctypedef struct Geometry:
        # Parameters part of the image geometry
        int   nVoxelX, nVoxelY, nVoxelZ;
        float sVoxelX, sVoxelY, sVoxelZ;
        float dVoxelX, dVoxelY, dVoxelZ;
        float *offOrigX;
        float *offOrigY;
        float *offOrigZ;
        float *DSO;

        # Parameters  of the Detector.
        int   nDetecU, nDetecV;
        float sDetecU, sDetecV;
        float dDetecU, dDetecV;
        float *offDetecU;
        float *offDetecV;
        float *DSD;
        float* dRoll;
        float* dPitch;
        float* dYaw;

        # The base unit we are working with in mm.
        float unitX;
        float unitY;
        float unitZ;

        #projection angle
        float alpha;
        float theta;
        float psi;
        # Centre of Rotation correction.
        float* COR;

        #Maximum length of cube
        float maxLength;

        #User option
        float accuracy;


    ctypedef struct Point3D:
        float x;
        float y;
        float z;

#TODO: Change from inline to stop duplication. Discussed on slack "python-questions" 28-Mar-2017
cdef inline void free_c_geometry(Geometry* c_geom):
    free(c_geom.offOrigX)
    free(c_geom.offOrigY)
    free(c_geom.offOrigZ)

    free(c_geom.offDetecU)
    free(c_geom.offDetecV)

    free(c_geom.dRoll)
    free(c_geom.dPitch)
    free(c_geom.dYaw)

    free(c_geom.COR)

    free(c_geom)

# Convert python-geometry to c-geometry.
# Interprets the size of np.array's in the C/cuda context.
# python-geometry's geometry.convert_contig_mode() is merged into this function. 
#TODO: Change from inline to stop duplication. Discussed on slack "python-questions" 28-Mar-2017
cdef inline Geometry* convert_to_c_geometry(p_geometry, int total_projections):
    cdef Geometry* c_geom =<Geometry *>malloc(sizeof(Geometry))

    ### Image ###
    c_geom.nVoxelX = p_geometry.nVoxel[2]
    c_geom.nVoxelY = p_geometry.nVoxel[1]
    c_geom.nVoxelZ = p_geometry.nVoxel[0]

    c_geom.sVoxelX = p_geometry.sVoxel[2]
    c_geom.sVoxelY = p_geometry.sVoxel[1]
    c_geom.sVoxelZ = p_geometry.sVoxel[0]
    
    c_geom.dVoxelX = p_geometry.dVoxel[2]
    c_geom.dVoxelY = p_geometry.dVoxel[1]
    c_geom.dVoxelZ = p_geometry.dVoxel[0]

    # TODO: array of constant for each alpha
    c_geom.offOrigX =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.offOrigX:
        raise MemoryError()
    c_geom.offOrigY =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.offOrigY:
        raise MemoryError()
    c_geom.offOrigZ =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.offOrigZ:
        raise MemoryError()
    c_geom.DSO =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.DSO:
        raise MemoryError()
    for i in range (total_projections):
        c_geom.offOrigX[i] = p_geometry.offOrigin[i][2]
        c_geom.offOrigY[i] = p_geometry.offOrigin[i][1]
        c_geom.offOrigZ[i] = p_geometry.offOrigin[i][0]
    for i in range(total_projections):
        c_geom.DSO[i] = p_geometry.DSO[i]

    ### Detector ###
    c_geom.nDetecU=p_geometry.nDetector[1]
    c_geom.nDetecV=p_geometry.nDetector[0]

    c_geom.sDetecU=p_geometry.sDetector[1]
    c_geom.sDetecV=p_geometry.sDetector[0]

    c_geom.dDetecU=p_geometry.dDetector[1]
    c_geom.dDetecV=p_geometry.dDetector[0]

    # TODO: array of constant for each alpha
    c_geom.offDetecU =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.offDetecU:
        raise MemoryError()
    c_geom.offDetecV =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.offDetecV:
        raise MemoryError()
    c_geom.DSD =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.DSD:
        raise MemoryError()
    for i in range (total_projections):
        c_geom.offDetecU[i] = p_geometry.offDetector[i][1]
        c_geom.offDetecV[i] = p_geometry.offDetector[i][0]
    for i in range(total_projections):
        c_geom.DSD[i] = p_geometry.DSD[i]

    # TODO: array of 0 for each alpha
    c_geom.dRoll =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.dRoll:
        raise MemoryError()
    c_geom.dPitch =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.dPitch:
        raise MemoryError()
    c_geom.dYaw =<float *>malloc(total_projections * sizeof(float))
    if not c_geom.dYaw:
        raise MemoryError()
    for i in range (total_projections):
        c_geom.dRoll[i] = p_geometry.rotDetector[i][2]
        c_geom.dPitch[i] = p_geometry.rotDetector[i][1]
        c_geom.dYaw[i] = p_geometry.rotDetector[i][0]

    # The base unit we are working with in mm.
    c_geom.unitX = 1
    c_geom.unitY = 1
    c_geom.unitZ = 1

    #projection angle
    # float alpha; #TODO: Check this is redundant

    # TODO: check this is correct line ~400 Ax.cpp. Seems like it can be omitted? Possibly redundant
    # Centre of Rotation correction.
    c_geom.COR =<float *>malloc(total_projections * sizeof(float))
    for i in range (total_projections):
        c_geom.COR[i] = p_geometry.COR[i]

    #Maximum length of cube
    # float maxLength; #TODO: Check this is redundant

    #User option
    c_geom.accuracy = p_geometry.accuracy

    return c_geom
