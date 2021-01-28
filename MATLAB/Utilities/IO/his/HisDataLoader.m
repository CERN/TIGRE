function [proj,geo,angles]=HisDataLoader(filepath)


if ~strcmp(filepath(end),'/') || ~strcmp(filepath(end),'\')
    filepath=[filepath,'/'];
end
%get all dicom files in this folder
x=dir([filepath, '*.his']);

proj(:,:,length(x))=readHISfile([filepath,x(end).name]);
for ii=1:length(x)-1
    proj(:,:,ii)=readHISfile([filepath,x(ii).name]);
end

end