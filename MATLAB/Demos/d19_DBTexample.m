%% DEMO 19: Example of DBT geometry
%
%   This is a simple example of using a DBT geometry in TIGRE. The DBT
%   geometry is based on the commercial GE equipment. We approximated the
%   currently CT geometry to the DBT. Please compare the results with other
%   toolboxes. Make sure to cite it properly.
%   You can compare this data with Lavi or FDA DBT toolboxes.
%
%   FDA -> https://github.com/DIDSR/ReconDBT
%   Lavi -> https://github.com/LAVI-USP/DBT-Reconstruction
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Rodrigo Vimieiro and Ander Biguri 
%--------------------------------------------------------------------------
%%  Geometry definition
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
%% Geometry structure:
%           -nVoxel:        3x1 matrix of number of voxels in the image
%           -sVoxel:        3x1 matrix with the total size in mm of the image
%           -dVoxel:        3x1 matrix with the size of each of the voxels in mm
%           -nDetector:     2x1 matrix of number of voxels in the detector plane
%           -sDetector:     2x1 matrix with the total size in mm of the detector
%           -dDetector:     2x1 matrix with the size of each of the pixels in the detector in mm
%           -DSD:           1x1 scalar value. Distance Source Detector, in mm
%           -DSO:           1x1 scalar value. Distance Source Origin.
%           -offOrigin:     3x1 matrix with the offset in mm of the centre of the image from the origin.
%           -offDetector:   2x1 matrix with the offset in mm of the centre
%           of the detector from the x axis
clc;clear;close all
%% Example
%
%
%
% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
% Distances
geo.DSD = 660;                              % Distance Source Detector      (mm)
geo.DSO = 620;                              % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[3062;2394];					% number of pixels              (px)
geo.dDetector=[0.1; 0.1]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[127;2053;784];                  % number of voxels              (vx)
geo.sVoxel=[63.3;205.3;78.4];               % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
Airgap = 22;                                % DBT airgap (mm)
geo.offOrigin =[((geo.sVoxel(1)/2)-...
(geo.DSD-geo.DSO)+Airgap);0;geo.sVoxel(3)/2]; % Offset of image from origin (mm)             
geo.offDetector=[0; geo.sDetector(2)/2];    % Offset of Detector            (mm)
                                            % These two can be also defined
                                            % per angle

% Auxiliary 
geo.accuracy=0.5;                           % Variable to define accuracy of
                                            % 'interpolated' projection
                                            % It defines the amoutn of
                                            % samples per voxel.
                                            % Recommended <=0.5             (vx/sample)

% Optional Parameters
% There is no need to define these unless you actually need them in your
% reconstruction
                                            
                                            
geo.COR=0;                                  % y direction displacement for 
                                            % centre of rotation
                                            % correction                   (mm)
                                            % This can also be defined per
                                            % angle

nprojs = 9;                                 % Number of projections
tubeangle = 25;                             % Angle range
angles=deg2rad(linspace(-tubeangle/2,tubeangle/2,nprojs)); % Angles in rad

geo.rotDetector=[0;0;0];                    % Rotation of the detector, by 
                                            % X,Y and Z axis respectively. (rad)
                                            % This can also be defined per
                                            % angle        
                                            
geo.mode='cone';                            % Or 'parallel'. Geometry type.

%% Plot the geometry                                           
plotgeometry(geo,0); 

%% Adapt CT geo to DBT
geo=staticDetectorGeo(geo,angles);

%% Adapt DBT projections to TIGRE CT projections
 
% Example of data (CT Head). This is not a true DBT, but it works as a
% example in TIGRE. Make sure to use a true tomosynthesis data.
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');

% If you use a true DBT projection, use the following lines to adapt you
% data to TIGRE CT.Remember to use -log in your data.
% projections = rot90(permute(projections,[2 1 3]),2);
% projections = -log(projections./single(2^14-1)); %(Use or not)

%% Simple BP
imgRecon = Atb( projections,geo,angles);
imgRecon = permute(imgRecon,[2 3 1]); 
imgRecon = rot90(imgRecon,2);

%% SART
imgRecon = SART(projections,geo,angles,2,'OrderStrategy','ordered');
imgRecon = permute(imgRecon,[2 3 1]);
imgRecon = rot90(imgRecon,2);

%%
function geo=staticDetectorGeo(geo,angles)
% This fucntion computes the translation and rotations needed in the system
% to desribe a detector that does not move when the source moves.

R=(geo.DSD-geo.DSO); %Radious of rotation of detector
geo.DSD=geo.DSD+(R-R*cos(angles));  % How much the detector is moving along X direction (Rodrigo comment)
geo.offDetector=[R*sin(angles); repmat(geo.offDetector(2),1,size(angles,2))]; % How much the detector is moving along Y direction (Rodrigo comment)
geo.rotDetector=[zeros(1,size(angles,2));zeros(1,size(angles,2));-angles]; % Rotates the detector in the opposite movement of the tube. (Rodrigo comment)

end
