function [geo, ScanXML,ReconXML] = ReadXimXML(Ximfolder)
% Dependence: xml2struct.m
% Date: 2020-03-28

%% Scan.xml
filestr = dir([Ximfolder filesep 'Scan.xml']);

ScanXML = getfield(xml2struct(fullfile(filestr.folder, filestr.name)),...
    'Scan');

%% Reconstruction.xml
filestr = dir([Ximfolder filesep '**' filesep 'Reconstruction.xml']);

ReconXML = getfield(xml2struct(fullfile(filestr.folder, filestr.name)),...
    'Reconstruction');

%% Generate geometry structure for TIGRE reconstruction
% Cone-beam CT scenario
geo.mode = 'cone';

%% ------------------- Load each parameter from ScanXML/ReconXML ---------------
% VARIABLE                                   DESCRIPTION                    UNITS
%-------------------------------------------------------------------------------------
% Distances
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[512; 512];					% number of pixels              (px)
geo.dDetector=[0.8; 0.8]; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[256;256;256];                   % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)              
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
                                            % These two can be also defined
                                            % per angle

% Auxiliary 
geo.accuracy=0.5;                           % Variable to define accuracy of
                                            % 'interpolated' projection
                                            % It defines the amoutn of
                                            % samples per voxel.
                                            % Recommended <=0.5             (vx/sample)
end