function [proj,geo,angles,dicomhdr]= dicomDataLoader(filepath)

% This file is a modification of code written for the following article:
%
% Hatamikia, S., Biguri, A., Kronreif, G., Kettenbach, J., Russ, T., Furtado, H.,
% Shiyam Sundar, L.K., Buschmann, M., Unger, E., Figl, M., Georg, D.
% and Birkfellner, W. (2020), 
% Optimization for customized trajectories in cone beam computed tomography. 
% Med. Phys., 47: 4786-4799. https://doi.org/10.1002/mp.14403
%
% Please cite it, especially if you are using C-arm CBCT

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

% Standard circular scans will only create 1 DICOM file, while one may
% "hack" the scan to do one-shots in each projection, and create several
% 1-projection scans. We do want to support both, so....:

if length(x)>1
    
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
else
    extra_required_fields={ 'PositionerPrimaryAngleIncrement',
                            'PositionerSecondaryAngleIncrement'};
    [dicomhdr,~]=checkDicomHeaders(dicomhdr,extra_required_fields);
    if isempty(dicomhdr)
        warning('DICOM header lacks required information to construct a geometry, outputing projection data only');
        dicomhdr{1}=dicominfo([filepath, x(1).name]);
        proj=dicomread(dicomhdr{1}.Filename);
        geo=[];
        angles=[];
        return;
    end
    % Load projetion data
    projaux=dicomread(dicomhdr{1}.Filename);
    if size(projaux,4)>1
        proj=single(squeeze(projaux(:,:,1,:)));
    else
        proj=single(projaux);
    end
    DSD=dicomhdr{1}.DistanceSourceToDetector;
    DSO=dicomhdr{1}.DistanceSourceToPatient;
    primary_angle=dicomhdr{1}.PositionerPrimaryAngle+dicomhdr{1}.PositionerPrimaryAngleIncrement'; %RAO/LAO
    secondary_angle=dicomhdr{1}.PositionerSecondaryAngle+dicomhdr{1}.PositionerSecondaryAngleIncrement'; %CRA/CAU
    % Check if projection has been padded.
    padding_fields={'ShutterLeftVerticalEdge',
                    'ShutterRightVerticalEdge',
                    'ShutterUpperHorizontalEdge',
                    'ShutterLowerHorizontalEdge'};
    if ~isempty(checkDicomHeaders(dicomhdr,padding_fields))
       % padding fields present, pad. 
       proj=proj(dicomhdr{1}.ShutterUpperHorizontalEdge:dicomhdr{1}.ShutterLowerHorizontalEdge,dicomhdr{1}.ShutterLeftVerticalEdge:dicomhdr{1}.ShutterRightVerticalEdge,:);
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
% Fill in the geometry info.
geo.dDetector=dicomhdr{1}.ImagerPixelSpacing;
geo.nDetector=double([size(proj,2);size(proj,1)]);
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
if all(secondary_angle==0)
    angles=-primary_angle*pi/180;
else
    angles=[-primary_angle;-secondary_angle;zeros(1,length(primary_angle))]*pi/180;
end

disp('Projection data loaded.')
end
% Function to check if dicomhdr has required fields
function [dicomhdr,correct_index]=checkDicomHeaders(dicomhdr,extra_required)

required_fields={'DistanceSourceToDetector',
    'DistanceSourceToPatient',
    'Width',
    'Height',
    'PositionerPrimaryAngle',
    'PositionerSecondaryAngle',
    'ImagerPixelSpacing'};
if nargin>1
    required_fields=cat(1,required_fields,extra_required);
end
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