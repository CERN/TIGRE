%% ALGOTESTS
% This function tests different algorithms of image recosntruction
% using the CUDA CBCT code.
% 
% Add new test of functionallity here, for future reference.


%% Set up enviroment.

clear;
clc;
close all;

% Initialize toolbox
initTOOLBOX;
%% Create geometry

% Load a default geometry (look inside if curious)
GeometrySettingRandoHead;
Geometry.nVoxel=[512;512;512]/4;
Geometry.sVoxel=[256;256;256];
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;

% Set accuracy (1/amount of samples per voxel)
Geometry.accuracy=1;

% Set projection angles
init=0;          % start angle
step=10;        % step
finish=360-step; % end angle

alpha=[init:step:finish]*pi/180;

%% Plot the geometry

plotgeometry(Geometry,30)

%% Define the test image.
 
% %%% 3D Shepp-Logan
img=phantom3dAniso('Modified Shepp-Logan',Geometry.nVoxel);

% %%% thorax digital phantom
% img=thoraxPhantom(Geometry.nVoxel);

%% Plot the image

plotImg(img,'Step',1,'Dim',3);

%% Generate some data

data=Ax(img,Geometry,alpha);
% Add noise to data;

maxval=max(data(:));
data_noise=imnoise(data/maxval,'poisson');
data_noise=data_noise.*maxval;
data=data_noise;
clear data_noise
%% Plot the data

plotProj(data,alpha);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niter=75;
%% SIRT 
[resSIRT, errSIRT]=SIRT_CBCT(data,Geometry,alpha,15);
%% SART

% Due to memory management stuff, SART is about 5 times slower
[resSART, errSART]=SART_CBCT(data,Geometry,alpha,15);

%% OS-SART
block_size=10;
[resOSSART, errOSSART]=OS_SART_CBCT(data,Geometry,alpha,10,'BlockSize',15);

%% FDK
[resFDK, errFDK]=FDK_CBCT(data,Geometry,alpha);

%% comparison

hold on
plot(errSIRT);
plot(errSART);
plot(errOSSART);
legend({'SIRT','SART','OS-SART'})
ylabel('L2 norm');
xlabel('iterations')
xlim([1 niter])

%%
figure('Name','Diferent reconstruction algorithms RANDO-HEAD, full data, 200 iterations')
slices=70;slices2=64;
ax1=subplot(141);


imgplot=squeeze(resSIRT(slices,:,:));
color=prctile(imgplot(:),[1 99]);
imshow([imgplot';resSIRT(:,:,slices2)],[],'Border','tight');axis xy; caxis(color);caxis([0, 0.35]);
title(['OS-SART '])


ax2=subplot(142);
imgplot=squeeze(resSART(slices,:,:));
color=prctile(imgplot(:),[1 99]);
imshow([imgplot';resSART(:,:,slices2)],[],'Border','tight');axis xy; caxis(color);caxis([0, 0.35]);
title(['OS-SART->init multigrid '])

ax3=subplot(143);
imgplot=squeeze(resOSSART(slices,:,:));
color=prctile(imgplot(:),[1 98]);
imshow([imgplot';resOSSART(:,:,slices2) ],[],'Border','tight');axis xy; caxis(color);caxis([0, 0.35]);
title(['OS-SART->init FDK '])

ax4=subplot(144);
imgplot=squeeze(resFDK(slices,:,:));
color=prctile(imgplot(:),[1 99]);
imshow([imgplot';resFDK(:,:,slices2) ],[],'Border','tight');axis xy; caxis(color);caxis([0, 0.35]);
title(['FDK '])
linkaxes([ax1,ax2,ax3,ax4], 'xy');
break;
