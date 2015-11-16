
% GEOMETRY: an structure describing the geometry of the CBCT, containing:
%           -nVoxel:        3x1              matrix of number of voxels in the image
%           -sVoxel:        3x1              matrix with the total size in mm of the image
%           -dVoxel:        3x1              matrix with the size of each of the voxels in mm
%           -nDetector:     2x1              matrix of number of voxels in the detector plane
%           -sDetector:     2x1              matrix with the total size in mm of the detector
%           -dDetector:     2x1              matrix with the size of each of the pixels in the detector in mm
%           -DSD:           1x1              scalar value. Distance Source Detector, in mm  
%           -DSO:           1x1              scalar value. Distance Source Origin.
%           -offOrigin:     3x1 or 3xNangles matrix with the offset in mm of the centre of the image from the origin.
%           -offDetector:   2x1 or 3xNangles matrix with the offset in mm of the centre of the detector from the x axis
%%  GEOMETRY DEFINITION
%                 
%                  Detector plane, behind
%              |-----------------------------| 
%              |                             | 
%              |                             | 
%              |                             | 
%  Centered    |                             | 
%    at O      A V    +--------+             |
%              |     /        /|             |
%     A Z      |    /        / |*D           |
%     |        |   +--------+  |             |
%     |        |   |        |  |             |
%     |        |   |     *O |  +             |
%     *--->y   |   |        | /              |
%    /         |   |        |/               |
%   V X        |   +--------+        U       |
%              .--------------------->-------|
%     
%            *S
%  
%  
%%
Geometry.DSD = 1536;   
Geometry.DSO = 1000;

Geometry.nDetector=[512; 512];
Geometry.dDetector=[0.8; 0.8];
Geometry.sDetector=Geometry.nDetector.*Geometry.dDetector;

Geometry.nVoxel=[256;256;256];
Geometry.sVoxel=Geometry.nVoxel; 
Geometry.dVoxel=[1; 1; 1];

Geometry.offOrigin=[0;0;0];           
Geometry.offDetector=[0; 0];