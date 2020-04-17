clear,clc
%% datasetfolder
datafolder = 'E:\BigData\Edge\CBCT_Export\2020-01-10_100318';
tic
%% Data loading: proj, angles geo
[proj, angles, geo] = EasyDataLoader(datafolder);

%% FDK recon:
angles = flip(angles);
if(abs(geo.offDetector(1)) > geo.sDetector(1)/3)
    warning("Half-fan FDK is applied\n");    
    imgFDK = FDK_half(proj, geo, angles);
else
    imgFDK = FDK(proj, geo, angles);
end
toc
%% Show results
imgFDK(imgFDK<0)=0;

%% plotImg(imgFDK./max(imgFDK(:)),'Dim','Z');
imshow(flipud(imgFDK(:,:,46)),[])

%%
imgFDK = ROIMask(imgFDK, 510);
imgFDK = imgFDK./max(imgFDK(:));

