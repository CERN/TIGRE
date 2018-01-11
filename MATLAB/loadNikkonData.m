function [proj,angles]=loadNikkonData(fpath,angles,nsamples,whitelevel)

%% compute index of desired data
maxdata=length(angles);
steps=floor(maxdata/nsamples);
idx=1:steps:steps*nsamples;

%% get filename
% assuming TIF and 4 digits.

mylist=cellstr(ls(fpath));
firstfile=find(~cellfun('isempty', strfind(mylist,'tif'))',1);
firstfile=mylist{firstfile};
filename=firstfile(1:end-8);

%% load images
proj=[];
for ii=1:nsamples
     proj(:,:,ii)=single(imread([fpath, '\' filename,num2str(idx(ii),'%04d'),'.tif']));
end
%% Beer lambert
proj=-log(proj/single(whitelevel));
angles=angles(idx);