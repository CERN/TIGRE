clear,clc
%% folder
datafolder = 'E:\BigData\Edge\CBCT_Export\2019-09-03_115522';

%% Data loading: proj, angles geo
[proj, angles, geo] = EasyDataLoader(datafolder);

%% FDK recon:
imgFDK=FDK(proj,geo,angles);

%% Show results
imgFDK(imgFDK<0)=0;

plotImg(imgFDK./max(imgFDK(:)),'Dim','Z');
