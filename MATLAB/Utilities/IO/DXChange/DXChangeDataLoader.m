function [proj,geo,angles]=DXChangeDataLoader(fname)
%DXChangeDataLoader(filepath) Loads DXChange datasets into TIGRE standard
%
%   DXChangeDataLoader(fname OPT,VAL,...) uses options and values. 
%      These are options in case you don't want to load the entire
%      dataset, but only particular sets of projections. 

%%
% make sure its h5 file
[~, ~, fExt] = fileparts(fname);
if ~strcmp(fExt,'.h5')
    error('Filename is not .h5 file')
end
proj=single(h5read(fname,"/exchange/data/"));
white=single(h5read(fname,"/exchange/data_white/"));
dark=single(h5read(fname,"/exchange/data_dark/"));
white(white<0)=0;
dark(dark<0)=0;
white=mean(white,3);
dark=mean(dark,3);
dark=zeros(size(dark),"single");
proj=(proj-dark)./(white-dark);
proj(proj<0)=0.01;
proj=-log(proj/(max(proj(:))));
proj=permute(proj,[2 1 3]);


angles=h5read(fname,"/exchange/theta").'*pi/180;
geo=defaultGeometry('mode','parallel','nDetector',[size(proj,2),size(proj,1)]);

end