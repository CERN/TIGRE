clc
close all;
%% datasetfolder
datafolder = 'D:\Edge\CBCT_Export\2019-09-04_133218';
%% Data loading: proj, angles geo
[proj, angles, geo] = EasyDataLoader(datafolder, 1);

%% FDK recon:
% half-fan FDK will be merged with FDK later
if(abs(geo.offDetector(1)) > geo.sDetector(1)/3)
    warning("Half-fan FDK is applied\n");    
    imgFDK = FDK_half(proj, geo, angles);
else
    imgFDK = FDK(proj, geo, angles);
end

%% Set all negatives to zero
imgFDK(imgFDK<0)=0;

%% imshow with a ROI mask
imgFDK = ROIMask(imgFDK, size(imgFDK,1));
figure;
subplot(1,3,1),imshow(imgFDK(:,:,30),[])
subplot(1,3,2),imshow(imgFDK(:,:,45),[])
subplot(1,3,3),imshow(imgFDK(:,:,80),[])
