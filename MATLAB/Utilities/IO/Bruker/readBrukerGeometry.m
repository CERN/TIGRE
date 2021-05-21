function [geo,angles]=readBrukerGeometry(folder_or_file)

% Developed by A. Biguri

if endsWith(folder_or_file,'.log')
    % is the log file itself
    fid=fopen(folder_or_file);
    isfolder=false;
else
    % is the folder where it lives
    isfolder=true;
    file = dir([folder_or_file,'/*.log']); %
    if isempty(file)
        error(['No .log file found in folder: ', folder_or_file]);
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
geo.nDetector=[str2double(xtekctText{2}(strcmp('Number of Columns', xtekctText{1})));
    str2double(xtekctText{2}(strcmp('Number of Rows', xtekctText{1})))];
% Size of pixels in the detector
geo.dDetector=[str2double(xtekctText{2}(strcmp('Camera Pixel Size (um)', xtekctText{1})))/1000;
    str2double(xtekctText{2}(strcmp('Camera Pixel Size (um)', xtekctText{1})))/1000* str2double(xtekctText{2}(strcmp('CameraXYRatio', xtekctText{1})))];
% Total size of the detector
geo.sDetector=geo.nDetector.*geo.dDetector;

%% Offset of the detector:
geo.offDetector=[0;0];

%% Image information
% Number of pixel in the detector

% Size of each pixel
geo.dVoxel=[str2double(xtekctText{2}(strcmp('Image Pixel Size (um)', xtekctText{1})))/1000;
    str2double(xtekctText{2}(strcmp('Image Pixel Size (um)', xtekctText{1})))/1000;
    str2double(xtekctText{2}(strcmp('Image Pixel Size (um)', xtekctText{1})))/1000 ];
% Size of the image in mm
geo.nVoxel=[geo.nDetector(1);geo.nDetector(1);geo.nDetector(2)];
geo.sVoxel=geo.nVoxel.*geo.dVoxel;

geo.offOrigin=[0;0;0];
%% Global geometry
geo.DSO=str2double(xtekctText{2}(strcmp('Object to Source (mm)', xtekctText{1})));
geo.DSD=str2double(xtekctText{2}(strcmp('Camera to Source (mm)', xtekctText{1})));


% Ignoring prior image information
mag=geo.DSD/geo.DSO;
geo.dVoxel=[geo.dDetector(1);geo.dDetector(1);geo.dDetector(2)]/mag;
geo.nVoxel=[geo.nDetector(1);geo.nDetector(1);geo.nDetector(2)];
geo.sVoxel=geo.nVoxel.*geo.dVoxel;


%% Detector offset
if endsWith(folder_or_file,'.log')
    folder=folder_or_file(1:end-4);
else
    folder=folder_or_file;
end
file = dir([folder,'/*.csv']); %
if ~isempty(file)
    geo.offDetector=csvread([file(1).folder, '/', file(1).name],5,1).';
end
%% whitelevel

geo.whitelevel=2^str2double(xtekctText{2}(strcmp('Depth (bits)', xtekctText{1})));

%%
%% angles

anlge_step=str2double(xtekctText{2}(strcmp('Rotation Step (deg)', xtekctText{1})));
initial_angle=0;
n_angles=str2double(xtekctText{2}(strcmp('Number of Files', xtekctText{1})))-1;
angles=(initial_angle:anlge_step:(initial_angle+360))*pi/180;
assert(size(angles,2)==n_angles,'Assertion failed: Inconsistent data detected. Number of projections and angle information do not match\n');

