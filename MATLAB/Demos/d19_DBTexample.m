clc;clear;close all
%%
%
% Code by Rodrigo Vimiero.
% https://github.com/rodrigovimieiro
%% Digital brest Thomosynthesis demo
%
% Few tricks to the geometry are needed to fix the detector in place.
% OTherwise its basic tomgoraphy.
%


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
geo.nDetector=[3064;2396]/2;					% number of pixels              (px)
geo.dDetector=[0.1; 0.1]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[128;2054;784]/2;                  % number of voxels              (vx)
geo.sVoxel=[63.3;205.3;78.4]/2;               % total size of the image       (mm)
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

%% Adapt CT geo to DBT
geo=staticDetectorGeo(geo,angles);
%% Adapt DBT projections to TIGRE CT projections

% Example of data (CT Head). This is not a true DBT, but it works as a
% example in TIGRE. Make sure to use a true tomosynthesis data.
head=headPhantom(geo.nVoxel);
projections=Ax(head,geo,angles,'interpolated');
% If you use a true DBT projection, use the following lines to adapt you
% data to TIGRE CT.Remember to use -log in your data.
% projections = -log(projections./single(2^14-1)); %your maximum white
% level.
%% Simple BP
recFDK = FDK( projections,geo,angles);

%% SART
recSART = SART(projections,geo,angles,2,'OrderStrategy','ordered');

%% plot
plotImg([recFDK recSART],'dim','x','clims',[0 1])

