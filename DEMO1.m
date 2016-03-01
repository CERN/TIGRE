%% How to use Ax with CUDA:
% 
%
% Attention: The code has been written for a NVIDIA Tesla k40 graphic card,
% and most likely need to be tunned if it is wanted to use with another a
% different one. 
%
%
% This code solves geometrically the A*x of a linear system of equations
% for the Cone Beam Computed Tomography geometry. In other words, it
% computes a projection with a given image. It can be used for both forward
% projection (simulating the measured 2D projections) and for iterative
% inverse solving of the image, where in many cases (as in SART algorithm)
% the A*x term has to be computed several times. The benefit is that there
% is no need of actually computing the matrix A.
%
% The code makes use of graphic cards hardware implemented texture memory 
% and the multy CPU arquitechture of modern GPUs to simultaneously compute
% the projection in each detector pixel.
%
%
% A single call to Ax will compute a single projection (as for 05/08/2015).
% Detector array is plane.
%
% INPUT ARGUMETS: 
%       projection=Ax(image, geometry)
%  
% IMAGE: a 3D image of the current solution.
% GEOMETRY: an structure describing the geometry of the CBCT, containing:
%           -nVoxel:        3x1 matrix of number of voxels in the image
%           -sVoxel:        3x1 matrix with the total size in mm of the image
%           -dVoxel:        3x1 matrix with the size of each of the voxels in mm
%           -nDetector:     2x1 matrix of number of voxels in the detector plane
%           -sDetector:     2x1 matrix with the total size in mm of the detector
%           -dDetector:     2x1 matrix with the size of each of the pixels in the detector in mm
%           -DSD:           1x1 scalar value. Distance Source Detector, in mm  
%           -DSO:           1x1 scalar value. Distance Source Origin.
%           -offOrigin:     3x1 matrix with the offset in mm of the centre of the image from the origin.
%           -offDetector:   2x1 matrix with the offset in mm of the centre of the detector from the x axis
% ALPHA: the angle in radians of the projection
% 

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
%  
%  
%  
%% Initialize toolbox
clear;
clc;
InitToolbox;

%% Geometry

% Image
Geometry.nVoxel=[128;128;128];                              % Voxel size
Geometry.sVoxel=[460;460;460];                              % Image size in mm
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;           % Voxel size in mm

% Detecotr
Geometry.nDetector=[128;128];                               % Detector size   
Geometry.sDetector=[ 1024;800];                             % Detector size in mm
Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;  % Size of each detector pixel

% distances
Geometry.DSD = 1500;    % Distance Source to Detector
Geometry.DSO = 1100;    % Distance Source to Origin (of XYZ, or mid image)

Geometry.offOrigin=[0; 0; 0];     % Rigid motion of image in mm. Can be 3x1 or 3xN        
Geometry.offDetector=[0;0];       % Offset of detector in mm. Can be 2x1 or 2XN
Geometry.accuracy=1;              % Smaller number->more accurate. Not recommended to be bigger than 1

alpha=[0:1:359]*pi/180;           % Anlges of projection
%% Use a digital phantom

phantom=thoraxPhantom(Geometry.nVoxel);

%% Generate data


data=Ax(phantom,Geometry,alpha); % this uses GPU
% plot the projections jumping 1 and save the result as a gif
plotProj(data,alpha,'Step',2,'Savegif','DEMO1.gif');
%% Reconstruct image

% FDK
[resFDK,errFDK]=FDK_CBCT(data,Geometry,alpha);

% OS-SART with multigrid initialization, 200 iterations and 20 simultaneous
% block size updates

niter=200;
nblock=20;
% [resOSSART,errOSSART]=OS_SART_CBCT(data,Geometry,alpha,niter,'BlockSize',nblock,'Init','multigrid');

%% Plot the result image

% plot XY slices
plotImg(resFDK,'Dim','Z')

% plot XZ slices
plotImg(resOSSART,'Dim','X')

