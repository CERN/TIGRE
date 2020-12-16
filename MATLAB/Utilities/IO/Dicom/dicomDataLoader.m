function [proj,geo,angles,dicomhdr]= dicomDataLoader(filepath)


if ~strcmp(filepath(end),'/') || ~strcmp(filepath(end),'\')
    filepath=[filepath,'/'];
end
%get all dicom files in this folder
x=dir([filepath, '*.dcm']);
% grab the last one (allows preallocation)
dicomhdr{length(x)}=dicominfo([filepath, x(length(x)).name]);
% If its not Philips, we don't know if its supported. 
if ~contains(dicomhdr{length(x)}.Manufacturer,'Philips')
   warning(['This loader has only been tested fo Philips scanners, which your DICOM data does not seem to be \n',...
             'If you have data not supported by this loader that you want supported, contact tigre.toolbox@gmail.com\n',...
             'and we can work together in improving this code to load your data.']);
   disp('Attempting to read the data anyway....');
end
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

end
function [dicomhdr,correct_index]=checkDicomHeaders(dicomhdr)

required_fields={'DistanceSourceToDetector',
                'DistanceSourceToPatient',
                 'Width',
                 'Height',
                 'PositionerPrimaryAngle',
                 'PositionerSecondaryAngle'};
invalid_index=false(1,length(dicomhdr));
correct_index=1:length(dicomhdr);
for ii=1:length(dicomhdr)
    if ~all(isfield(dicomhdr{ii},required_fields))
        invalid_index(ii)=true;
    end
end
correct_index(invalid_index)=[];
dicomhdr(invalid_index)=[];
end