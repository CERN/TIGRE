%clear,clc
close all;

%% datasetfolder
datafolder = 'E:\BigData\Edge\CBCT_Export\2019-09-03_115522';
tic
%% Data loading: proj, angles geo
[proj, angles, geo] = EasyDataLoader(datafolder, 1);

%% FDK recon:
% img matrix rotation problem: hard-coded somewhere in FDK
if(abs(geo.offDetector(1)) > geo.sDetector(1)/3)
    warning("Half-fan FDK is applied\n");    
    imgFDK = FDK_half(proj, geo, angles);
else
    imgFDK = FDK(proj, geo, angles);
end
toc

if(angles(end)- angles(1)>0)
    runtitle = "Clockwise";
    fprintf("%s\n", runtitle);
else
    runtitle = "Counter-closewise";
    fprintf("%s\n", runtitle);
end
%% Show results
imgFDK(imgFDK<0)=0;

%% plotImg(imgFDK./max(imgFDK(:)),'Dim','Z');
imgFDK = ROIMask(imgFDK, 510);
figure;
subplot(1,3,1),imshow(imgFDK(:,:,30),[])
subplot(1,3,2),imshow(imgFDK(:,:,45),[])
subplot(1,3,3),imshow(imgFDK(:,:,80),[])
title(runtitle)
sound(rand(200,1))
