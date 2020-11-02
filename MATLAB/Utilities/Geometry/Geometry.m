classdef Geometry < GeometryInterface
% Describing your geometry
%
%  In TIGRE the geometry is stored in an structure. To see documentation
%  about geometry, run:
%     
%     doc('TIGRE/Geometry')
%
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
% Coded by:           Ander Biguri 
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

% Shows Geometry diagram
    properties
        % VARIABLE                                   DESCRIPTION                    UNITS
        %-------------------------------------------------------------------------------------
        % Distances
        DSD                                     % Distance Source Detector      (mm)
        DSO                                     % Distance Source Origin        (mm)
        % Detector parameters
        nDetector            					% number of pixels              (px)
        dDetector           					% size of each pixel            (mm)
        sDetector                               % total size of the detector    (mm)
        % Image parameters
        nVoxel                                  % number of voxels              (vx)
        sVoxel                                  % total size of the image       (mm)
        dVoxel                                  % size of each voxel            (mm)
        % Offsets
        offOrigin                               % Offset of image from origin   (mm)              
        offDetector                             % Offset of Detector            (mm)
                                                % These two can be also defined
                                                % per angle
                                                
        % Auxiliary 
        accuracy=0.5;                           % Variable to define accuracy of
                                                % 'interpolated' projection
                                                % It defines the amoutn of
                                                % samples per voxel.
                                                % Recommended <=0.5             (vx/sample)

        % Optional Parameters
        % There is no need to define these unless you actually need them in your
        % reconstruction


        COR=0;                                  % y direction displacement for 
                                                % centre of rotation
                                                % correction                   (mm)
                                                % This can also be defined per
                                                % angle

        rotDetector=[0;0;0];                    % Rotation of the detector, by 
                                                % X,Y and Z axis respectively. (rad)
                                                % This can also be defined per
                                                % angle        

        mode='cone';                            % Or 'parallel'. Geometry type. 
    end
end
