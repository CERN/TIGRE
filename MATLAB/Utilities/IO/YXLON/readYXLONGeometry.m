function [geo,angles]=readYXLONGeometry(folder)

% Developed by A. Biguri and W. Sun
% W. Sun edited on 06.10.2018 for 3D pro version 5.2.6809.15380 (2018)
% Modified for YXLON by A. Biguri

file = dir([folder, '\ScanParameter.xml']);
if ~isempty(file)
    [geo,angles]=loadXMLGeometry(folder);
else
    
    file = dir([folder,'\*.txt']); %
    if isempty(file)
        error(['No .txt file found in folder: ', folder]);
    end
    [geo,angles]=loadTXTGeometry(folder);
end
end

function [geo,angles]=loadXMLGeometry(folder)

    s = readstruct([folder, '\ScanParameter.xml']);
    geodata=s.ConeBeamScanParameter;
    
    if isfield(geodata,'DetectorOrientation') && ~strcmpi(geodata.DetectorOrientation,'original')
       warning('Untested Detector orientation found') 
    end
    if isfield(geodata,'DetectorFlipDirections') && ~strcmpi(geodata.DetectorFlipDirections,'none')
       warning('Untested Detector flip direction found') 
    end
    
    geo.dDetector=[geodata.PixelSize.Width; geodata.PixelSize.Height];
    geo.nDetector=[geodata.ProjectionSize.Width; geodata.ProjectionSize.Height];
    geo.sDetector = geo.nDetector .* geo.dDetector;
    
    if isfield(geodata,'ScanTrajectory')
        if ~strcmpi(geodata.ScanTrajectory.xsi_typeAttribute,'CircleTrajectory')
            error(['Unsuported trajectory ' geodata.ScanTrajectory.xsi_typeAttribute])
        end
    else
       error('Unable to identify scan trajectory')
    end
    
    geo.mode='cone';
    scale=geodata.ScanTrajectory.PixelPitch;
    geo.DSD=geodata.ScanTrajectory.Fdd*scale;
    geo.DSO=geodata.ScanTrajectory.Fod*scale;
    geo.COR=-scale*geodata.ScanTrajectory.YoCenter;
    
    if isfield(geodata,'DejustmentParameter')
        geo.offDetector=-scale*[geodata.DejustmentParameter.HorizontalDetectorOffset; geodata.DejustmentParameter.VerticalDetectorOffset];

        if geodata.DejustmentParameter.DetectorTiltA~=0 || geodata.DejustmentParameter.DetectorTiltB~=0 || geodata.DejustmentParameter.DetectorTiltC~=0
            warning("Detector tilt detected but not accounted for. If you have more info on this email tigre.toolbox@gmail.com")
        end
    end
    
    geo.nVoxel=[geo.nDetector(1);geo.nDetector(1);geo.nDetector(2)];
    geo.sVoxel=[geo.sDetector(1);geo.sDetector(1);geo.sDetector(2)]*geo.DSO/geo.DSD;
    geo.dVoxel=geo.sVoxel./geo.nVoxel;
    geo.whitelevel=2^geodata.ProjectionBitdepth-1;
    
    initial_angle=geodata.ScanTrajectory.StartAngle;
    n_angles=geodata.ScanTrajectory.NumberOfProjections;
    range_angles=geodata.ScanTrajectory.Rotation*2*pi/360;
    angles=linspace(initial_angle,initial_angle+range_angles,n_angles+1);
    angles(end)=[];
    
end

function [geo,angles]=loadTXTGeometry(folder)

    file = dir([folder,'\*.txt']); %
    if isempty(file)
        error(['No .txt file found in folder: ', folder]);
    end
    filename=file(1).name;
    fid=fopen([folder,'/',filename]);
% is the folder where it lives
file = dir([folder_or_file,'\*.txt']); %
if isempty(file)
    error(['No .txt file found in folder: ', folder_or_file]);
end
filename=file(1).name;
fid=fopen([folder_or_file,'/',filename]);

% check if it was opened correctly
if(fid==-1)
    error('Wrong file path');
end
ETCText = textscan(fid, '%q %q','Delimiter','\t','MultipleDelimsAsOne',1);
fclose(fid);

%% Detector information
%just to create a BASIC geometry without garbage
% Number of pixels in the detector
geo.nDetector=[str2double(ETCText{1}(strcmp('Number of Detector Elements ', ETCText{2})));
    str2double(ETCText{1}(strcmp('Number of Detector Elements ', ETCText{2})))];

% Size of pixels in the detector
geo.dDetector=[str2double(ETCText{1}(strcmp('2D-Pixel Size [mm]', ETCText{2})));
    str2double(ETCText{1}(strcmp('2D-Pixel Size [mm]', ETCText{2})))];

geo.dDetector=geo.dDetector.*str2double(ETCText{1}(strcmp('Pixelbinning', ETCText{2})));
% Total size of the detector
geo.sDetector=geo.nDetector.*geo.dDetector;

%% Image information
% geo.nVoxel = [632;640;640];
geo.nVoxel=[str2double(ETCText{1}(strcmp('Image Dimension', ETCText{2})));
    str2double(ETCText{1}(strcmp('Image Dimension', ETCText{2})));
    str2double(ETCText{1}(strcmp('Number of Z-slices', ETCText{2}))) ];

% Size of each pixel
geo.dVoxel=[str2double(ETCText{1}(strcmp('3D-XY-Pixel Size [mm]', ETCText{2})));
    str2double(ETCText{1}(strcmp('3D-XY-Pixel Size [mm]', ETCText{2})));
    str2double(ETCText{1}(strcmp('3D-Z-Pixel Size [mm]', ETCText{2}))) ];

% Size of the image in mm
geo.sVoxel=geo.nVoxel.*geo.dVoxel;

%SDD and SOD
geo.DSO=str2double(ETCText{1}(strcmp('FOD [mm]', ETCText{2})));
geo.DSD=str2double(ETCText{1}(strcmp('FDD [mm]', ETCText{2})));

%extra bits
geo.offDetector=[0;0];
geo.offOrigin=[0;0;0];
geo.COR=str2double(ETCText{1}(strcmp('Center Offset [mm]', ETCText{2})));;

geo.accuracy=0.5;
geo.mode='cone';


geo.whitelevel=2^16-1;

%% angles
% 1st attempt
angles=linspace(0,2*pi,(str2double(ETCText{1}(strcmp('CT: Number of Projections', ETCText{2}))))+1);

angles=angles(1:length(angles)-1);



end