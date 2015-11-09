// Assume non-normalized geometry
__global__ void kernelPixelDetector( Geometry geo,
                                    double* detector,Point3D S){
   
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nDetecU* geo.nDetecV)
        return;
    
    
   

    /////////////////////////////// Get pixel coords
    int pixelU = idx % geo.nDetecU;
    int pixelV = idx / geo.nDetecU;
    
    //Start point
    Point3D S;
    S.x=-geo.DSO;
    S.y=0;
    S.z=0;
    
    //End point
    Point3D P;
    P.x=-(geo.DSD-geo.DSO);
    P.y= geo.dDetecU*(pixelU-(double)(geo.nDetecU/2)+0.5);
    P.z= geo.dDetecV*((double)(geo.nDetecV/2)+0.5-pixelV);
    
    // Geomtric trasnformations:
    
    //1: Offset detector
       
    //P.x
    P.y=P.y+geo.offDetectU;
    P.z=P.z+geo.offDetectV;
    //S doesnt need to chagne
    
    //2: Offset image (instead of offseting image, -offset everything else)
    
    P.x=P.x-geo.offOrigX;
    P.y=P.y-geo.offOrigY;
    P.z=P.z-geo.offOrigZ;
    
    S.x=S.x-geo.offOrigX;
    S.y=S.y-geo.offOrigY;
    S.z=S.z-geo.offOrigZ;
    
    //3: Rotate (around z)!
    Point3D P2;   
    P2.x=P.x*cos(geo.aplha)-P.y*sin(geo.alpha);
    P2.y=P.y*cos(geo.aplha)+P.x*sin(geo.alpha);
    P2.z=P.z;
    Point3D S2; 
    S2.x=S.x*cos(geo.aplha)-S.y*sin(geo.alpha);
    S2.y=S.y*cos(geo.aplha)+S.x*sin(geo.alpha);
    S2.z=S.z;
    // As we want the (0,0,0) to be in a corner of the image, we need to translate everything (after rotation);
    P2.x=P2.x-geo.sVoxelX/2;
    P2.y=P2.y-geo.sVoxelY/2;
    P2.z=P2.z-geo.sVoxelZ/2;
    S2.x=S2.x-geo.sVoxelX/2;
    S2.y=S2.y-geo.sVoxelY/2;
    S2.z=S2.z-geo.sVoxelZ/2;
    
    //4. Scale everything so dVoxel==1
    P2.x=P2.x/geo.dVoxelX;
    P2.y=P2.y/geo.dVoxelY;
    P2.z=P2.z/geo.dVoxelZ;
    S2.x=S2.x/geo.dVoxelX;
    S2.y=S2.y/geo.dVoxelY;
    S2.z=S2.z/geo.dVoxelZ;
    
   
    

    double vectX,vectY,vectZ;
    vectX=(P2.x-S2.x)/(geo.maxLength); 
    vectY=(P2.y-S2.y)/(geo.maxLength); 
    vectZ=(P2.z-S2.z)/(geo.maxLength);
    //here comes the deal
    double x,y,z;
    double sum=0;
    double i;

    for (i=0; i<=geo.maxLength; i=i+1){
        x=vectX*(double)i+S2.x;
        y=vectY*(double)i+S2.y;
        z=vectZ*(double)i+S2.z;
        // Make sure we have an image for the index we are looking for. Check if out of image.
        //if(x<0 | y<0 | z<0 | x> (double)geo.nVoxelX-1.0 | y> (double)geo.nVoxelY-1.0 | z> (double)geo.nVoxelZ-1.0){
        //    continue;   
        //}
        sum += (double)tex3D(tex, x+0.5, y+0.5, z+0.5);
    }
    detector[idx]=sum*(sqrt((S.x-P.x)*(S.x-P.x)+(S.y-P.y)*(S.y-P.y)+(S.z-P.z)*(S.z-P.z))/geo.maxLength);
}





//Normalized geometry input
__global__ void kernelPixelDetector( Geometry geo,
                                    double* detector,Point3D S){
   
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx>= geo.nDetecU* geo.nDetecV)
        return;
    
    
   

    /////////////////////////////// Get pixel coords
    int pixelU = idx % geo.nDetecU;
    int pixelV = idx / geo.nDetecU;
    
 
    // Get XYZ pixel coords
    
    Point3D P;    
    P.x=-(geo.DSD-geo.DSO)/*DDO*/                                              *cos(geo.alpha)   -   
         (geo.dDetecU*(pixelU-(double)(geo.nDetecU/2)+0.5) + geo.offDetecU)    *sin(geo.alpha);
    
    P.y= (geo.dDetecU*(pixelU-(double)(geo.nDetecU/2)+0.5) + geo.offDetecU)    *cos(geo.alpha)  -
         (geo.DSD-geo.DSO)                                                     *sin(geo.alpha);    
    
    P.z=  geo.dDetecV*((double)(geo.nDetecV/2)+0.5-pixelV);
  
    //Trasnlate coordinates to a corner of the image.
    //, so Axes are aligned with image boudnaries and 0,0,0 is a corner of the image 
    //( there is anly a single way of acieving this!) 
    P.x=P.x+(double)geo.nVoxelX/(double)2;
    P.y=P.y+(double)geo.nVoxelY/(double)2;
    P.z=P.z+(double)geo.nVoxelZ/(double)2;
    
    

    double vectX,vectY,vectZ;
    vectX=(P.x-S.x)/(geo.maxLength); 
    vectY=(P.y-S.y)/(geo.maxLength); 
    vectZ=(P.z-S.z)/(geo.maxLength);
    //here comes the deal
    double x,y,z;
    double sum=0;
    double i;

    for (i=0; i<=geo.maxLength; i=i+1){
        x=vectX*(double)i+S.x;
        y=vectY*(double)i+S.y;
        z=vectZ*(double)i+S.z;
        // Make sure we have an image for the index we are looking for. Check if out of image.
        //if(x<0 | y<0 | z<0 | x> (double)geo.nVoxelX-1.0 | y> (double)geo.nVoxelY-1.0 | z> (double)geo.nVoxelZ-1.0){
        //    continue;   
        //}
        sum += (double)tex3D(tex, x+0.5, y+0.5, z+0.5);
    }
    detector[idx]=sum*(sqrt((S.x-P.x)*(S.x-P.x)+(S.y-P.y)*(S.y-P.y)+(S.z-P.z)*(S.z-P.z))/geo.maxLength);
}
