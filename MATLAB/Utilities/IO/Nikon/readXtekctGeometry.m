function [geo,angles]=readXtekctGeometry(folder_or_file)

% Developed by A. Biguri and W. Sun
% W. Sun edited on 06.10.2018 for 3D pro version 5.2.6809.15380 (2018)

if endsWith(folder_or_file,'.xtekct')
    % is the xtekct file itself
    fid=fopen(folder_or_file);
    isfolder=false;
else
    % is the folder where it lives
    isfolder=true;
    file = dir([folder_or_file,'/*.xtekct']); %
    if isempty(file)
        error(['No .xtekct file found in folder: ', folder_or_file]);
    end
    filename=file(1).name;
    fid=fopen([folder_or_file,'/',filename]);
end
% check if it was oppened right
if(fid==-1)
    error('Wrong file path');
end
xtekctText = textscan(fid, '%s %s', 'Delimiter', '=', 'HeaderLines', 1, 'CommentStyle', '[');
fclose(fid);

%% Detector information
% Number of pixel in the detector
geo.nDetector=[str2double(xtekctText{2}(strcmp('DetectorPixelsX', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('DetectorPixelsY', xtekctText{1})))];
% Size of pixels in the detector
geo.dDetector=[str2double(xtekctText{2}(strcmp('DetectorPixelSizeX', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('DetectorPixelSizeY', xtekctText{1})))];
% Total size of the detector
geo.sDetector=geo.nDetector.*geo.dDetector;

%% Offset of the detector:
geo.offDetector=[str2double(xtekctText{2}(strcmp('DetectorOffsetX', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('DetectorOffsetY', xtekctText{1})))  ];

%% Image information
% Number of pixel in the detector
geo.nVoxel=[str2double(xtekctText{2}(strcmp('VoxelsX', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('VoxelsY', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('VoxelsZ', xtekctText{1}))) ];
% Size of each pixel
geo.dVoxel=[str2double(xtekctText{2}(strcmp('VoxelSizeX', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('VoxelSizeY', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('VoxelSizeZ', xtekctText{1}))) ];
% Size of the image in mm
geo.sVoxel=geo.nVoxel.*geo.dVoxel;
geo.offOrigin=[0;0;0];
%% Global geometry
geo.DSO=str2double(xtekctText{2}(strcmp('SrcToObject', xtekctText{1})));
geo.DSD=str2double(xtekctText{2}(strcmp('SrcToDetector', xtekctText{1})));
geo.COR=-str2double(xtekctText{2}(strcmp('CentreOfRotationTop', xtekctText{1})));
if (geo.COR==0)
    warning('Centre of Rotation seems to be zero. Make sure that it is true and that the machine did not omit that information.');
else
    warning('TIGRE doesnt know if the sign of COR is the right one. Consider triying both and reporting to tigre.toolbox@gmail.com.');
end
%% whitelevel

geo.whitelevel=str2double(xtekctText{2}(strcmp('WhiteLevel', xtekctText{1})));

%%
%% angles
% It can be either an .ang  or .txt file
% .ang
if ~isfolder
    % rewrite it. 
    [folder_or_file,~,~] = fileparts(folder_or_file);
end

    
file_angles=dir([folder_or_file,'/','*.ang']);
if ~isempty(file_angles)
    fid=fopen([folder_or_file,'/',file_angles.name]);
    xtekctText = textscan(fid, '%s %s', 'Delimiter', '\t', 'HeaderLines', 1);
    fclose(fid);
    angles=str2double(xtekctText{2})*pi/180;
    angles=angles.';
    return
end
% ._ctdata.txt
file_angles=dir([folder_or_file,'/','_ctdata.txt']);
if ~isempty(file_angles)
    fid=fopen([folder_or_file,'/',file_angles.name]);
    xtekctText = textscan(fid, '%s %s %s', 'Delimiter', '\t', 'HeaderLines', 3);
    fclose(fid);
    angles=str2double(xtekctText{2})*pi/180;
    angles=angles.';
    return
end

warning('File with definition of angles not found, estimating them from geometri info.')
angle_step=str2double(xtekctText{2}(strcmp('AngularStep', xtekctText{1})));
initial_angle=str2double(xtekctText{2}(strcmp('InitialAngle', xtekctText{1})));
n_angles=str2double(xtekctText{2}(strcmp('Projections', xtekctText{1})));
angles=initial_angle:angle_step*pi/180:(initial_angle+(n_angles-1)*angle_step)*pi/180;
assert(size(angles,2)==n_angles,'Assertion failed: Inconsistent data detected. Number of projections and angle information do not match\n');

