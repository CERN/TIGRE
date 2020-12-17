function [proj,geo,angles,dicomhdr]= dicomDataLoader(filepath)


if ~strcmp(filepath(end),'/') || ~strcmp(filepath(end),'\')
    filepath=[filepath,'/'];
end
%get all dicom files in this folder
x=dir([filepath, '*.dcm']);
% grab the last one (allows preallocation)
dicomhdr{length(x)}=dicominfo([filepath, x(length(x)).name]);
% If its not Philips, we don't know if its supported.
if ~isfield(dicomhdr{length(x)},'Manufacturer') || ~contains(dicomhdr{length(x)}.Manufacturer,'Philips')
    warning(sprintf(['This loader has only been tested fo Philips scanners, which your DICOM data does not seem to be \n',...
        'If you have data not supported by this loader that you want supported, contact tigre.toolbox@gmail.com\n',...
        'and we can work together in improving this code to load your data.']));
    disp('Attempting to read the data anyway....');
end
if ~isfield(dicomhdr{length(x)},'PositionerType') || ~strcmi(dicomhdr{length(x)}.PositionerType,'CARM')
    warning(sprintf(['This loader has only been tested with C-arm (CARM) scanners, which your DICOM data does not seem to be\n'...
        'If you have data not supported by this loader that you want supported, contact tigre.toolbox@gmail.com\n',...
        'and we can work together in improving this code to load your data.']))
    disp('Attempting to read the data anyway....');
    
end
% We have only teste CARM, but
% load info
for ii=1:length(x)-1
    dicomhdr{ii}=dicominfo([filepath, x(ii).name]);
end

% Its not usual, but there is a chance that some of the headers are
% corrupted. We should check if the required info is there, and otherwise
% error.
[dicomhdr,~]=checkDicomHeaders(dicomhdr);

% Load projetion data
proj=zeros(dicomhdr{1}.Height,dicomhdr{1}.Width,length(dicomhdr),'single');
for ii=1:length(dicomhdr)
    projaux=dicomread(dicomhdr{ii}.Filename);
    if size(projaux,4)>1
        proj(:,:,ii)=projaux(:,:,1,1);
    else
        proj(:,:,ii)=projaux;
    end
end
if any(proj(:)==0)
    warning(sprintf(['Zeroes detected in projection data, adding +1 to apply the Beer-Lambert law\n',...
        'This is caise by either very high attenuation in an image, to the point of blocking x-rays totally or the detector being cropped digitally (horizontal/vertical zero bands)\n',...
        'In the first instance, you may consider scanning your data again, in the second instance you may consider cropping the projections and changing the geo.xDetector information.\n',...
        'Finally, this can be caused by simply a rogue pixel, and you may continue safely.']))
    proj=proj+1;
end
proj=-log(proj./(max(proj(:))+1));

% fill data
DSD=zeros(1,length(dicomhdr));
DSO=zeros(1,length(dicomhdr));
primary_angle=zeros(1,length(dicomhdr));
secondary_angle=zeros(1,length(dicomhdr));

for ii=1:length(dicomhdr)
    DSD(ii)=dicomhdr{ii}.DistanceSourceToDetector;
    DSO(ii)=dicomhdr{ii}.DistanceSourceToPatient;
    primary_angle(ii)=dicomhdr{ii}.PositionerPrimaryAngle; %RAO/LAO
    secondary_angle(ii)=dicomhdr{ii}.PositionerSecondaryAngle; %CRA/CAU
end

% Fill in the geometry info.
geo.dDetector=dicomhdr{1}.ImagerPixelSpacing;
geo.nDetector=double([dicomhdr{1}.Height;dicomhdr{1}.Width]);
geo.sDetector=geo.dDetector.*geo.nDetector;
geo.offDetector=[0;0];
if length(unique(DSD))==1
    geo.DSD=DSD(1);
else % very unlikely in a medical CT scan
    geo.DSD=DSD;
end
if length(unique(DSO))==1
    geo.DSO=DSO(1);
else % very unlikely in a medical CT scan
    geo.DSO=DSO;
end
disp('Estimating acceptable image size...')
geo.dVoxel=[geo.dDetector(1);geo.dDetector(1); geo.dDetector(2)]*geo.DSO/geo.DSD;
geo.nVoxel=ceil([geo.nDetector(1)+abs(geo.offDetector(1))/geo.dDetector(1);geo.nDetector(1)+abs(geo.offDetector(1))/geo.dDetector(1);geo.nDetector(2)]);
geo.sVoxel=geo.nVoxel.*geo.dVoxel;

% this is not strictly necesary, but will make users understand the
% geometry easier. 
% TODO: make PositionerSecondaryAngle not required at all... Need more
% data. to test
if all(secondary_angle(ii)==0)
    angles=primary_angle;
else
    angles=[-primary_angle;-secondary_angle;zeros(1,length(primary_angle))]*pi/180;
end
disp('Projection data loaded.')
end
% Function to check if dicomhdr has required fields
function [dicomhdr,correct_index]=checkDicomHeaders(dicomhdr)

required_fields={'DistanceSourceToDetector',
    'DistanceSourceToPatient',
    'Width',
    'Height',
    'PositionerPrimaryAngle',
    'PositionerSecondaryAngle',
    'ImagerPixelSpacing'};
invalid_index=false(1,length(dicomhdr));
correct_index=1:length(dicomhdr);
for ii=1:length(dicomhdr)
    if ~all(isfield(dicomhdr{ii},required_fields))
        invalid_index(ii)=true;
    end
end
correct_index(invalid_index)=[];
dicomhdr(invalid_index)=[];

if sum(invalid_index)
    warning(['Some dicom files did not contain required fields: ', num2str(sum(invalid_index)),' files ignored.'])
end
end