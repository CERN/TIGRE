%% test real data
%% Init
clear;
clc;
close all;

% Initialize toolbox
initTOOLBOX;

%% Set up geometry
Geometry.DSD = 1536;   
Geometry.DSO = 1000;

Geometry.nDetector=[512; 512];
Geometry.dDetector=[0.8; 0.8];
Geometry.sDetector=Geometry.nDetector.*Geometry.dDetector;

Geometry.nVoxel=[512;512;512];
Geometry.sVoxel=[256;256;256];
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;

Geometry.offOrigin=[0;0;0];           
Geometry.accuracy=1;

addpath('C:\VOL_CT_modified\')
[P,D] = xread('C:\VOL_CT_modified\rando_head\');

alpha=P.Angle';

 % TESTED, the offsets should be like this
Geometry.offDetector=[-P.Uoff'*Geometry.dDetector(1);-P.Voff'*Geometry.dDetector(2)];


%% load data
data=zeros([Geometry.nDetector' length(alpha)]);


for ii=1:size(D,3)
    RealProj=double(D{ii});
    I0 = max(RealProj(:)); % if you don't know I0 
    Proj = -log(RealProj/I0);
    data(:,:,ii)=flipud(Proj);
    data(end,:,ii)=data(end-1,:,ii);
end

% data=data(:,:,1:10:end); 
% alpha=alpha(1:10:end);
% clear D I0 RealProj Proj


%% visualize projections.
vis=0;
if vis
    plotProj(data,alpha);
end
close all


%% reconstruct OS-SART

niter=20;
Geometry.nVoxel=[512;512;512]/2;
Geometry.sVoxel=[256;256;256];
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;

% [resSART,errSART]=SART_CBCT(data,Geometry,alpha,niter);
% [resOSSART,errOSSART]=OS_SART_CBCT(data,Geometry,alpha,niter,20);
% [resSIRT,errSIRT]=SIRT_CBCT(data,Geometry,alpha,niter);
[resFDK,errFDK]=FDK_CBCT(data,Geometry,alpha);

[resmulti,multiL2]=OS_SART_CBCT(data,Geometry,alpha,niter,'BlockSize',20,'Init','multigrid');
[resnoinit,noinitL2]=OS_SART_CBCT(data,Geometry,alpha,niter,'BlockSize',20);
[resinitFDK,FDKL2]=OS_SART_CBCT(data,Geometry,alpha,niter,'BlockSize',20,'Init','FDK');
[resFDK,errFDK]=FDK_CBCT(data,Geometry,alpha);

save rando_head20_256_different_init.mat resmulti multiL2 resnoinit noinitL2 resinitFDK FDKL2 resFDK errFDK
break;
%% Plot errors
figure('Name','L2 errors')
% plot([errSART]);
hold on;
plot([errOSSART]);
% plot([errSIRT]);
plot([1,niter],[errFDK,errFDK])
legend({'SART','OS-SART','SIRT','FDK'})
%% Plot images

plotthis=0; % this is so you dont accidentally spend hours looking at the plots
if plotthis
    plotImg(resSART,'Step',2,'Dim',3);
    plotImg(resOSSART,'Step',2,'Dim',3);
    plotImg(resSIRT,'Step',2,'Dim',3);
    plotImg(resFDK,'Step',1,'Dim',3);
end




%% Just OS-SART and FDK

figure('Name','Diferent reconstruction algorithms RANDO-HEAD')
slices=256;



ax2=subplot(121);
imgplot=squeeze(resOSSART(slices,:,:));
color=prctile(imgplot(:),[1 99]);

imshow([imgplot';resOSSART(:,:,140)],[],'Border','tight');axis xy; caxis(color);
title(['OS-SART ', num2str(niter), ' iterations'])

ax1=subplot(122);
color=prctile(imgplot(:),[1 99]);
imgplot=squeeze(resFDK(slices,:,:));
color=prctile(imgplot(:),[1 99]);

imshow([imgplot';resFDK(:,:,140) ],[],'Border','tight');axis xy; caxis(color);
title(['FDK '])



linkaxes([ax1,ax2], 'xy');

break;

%% Compare images
figure('Name','Diferent reconstruction algorithms RANDO-HEAD')

subplot(221);

imgplot=squeeze(resFDK(slices,:,:));

color=prctile(imgplot(:),[1 99]);
imshow(imgplot',[]);axis xy; caxis(color);
title(['FDK '])



subplot(222);
slices=124;
imgplot=squeeze(resOSART(slices,:,:));
imshow(imgplot',[]);axis xy; caxis(color);
title(['SART: ', num2str(niter), ' iterations'])

subplot(223);

imgplot=squeeze(resOSSART(slices,:,:));
imshow(imgplot',[]);axis xy; caxis(color);
title(['OS-SART ', num2str(niter), ' iterations'])


subplot(224);

imgplot=squeeze(resSIRT(slices,:,:));
imshow(imgplot',[]);axis xy; caxis(color);
title(['SIRT ', num2str(niter), ' iterations'])



%% Random crap that I migth want to keep
%%
%%
%%
%% Post process
[x, y]=meshgrid(1:Geometry.nVoxel(1),1:Geometry.nVoxel(2));
x=x-(Geometry.nVoxel(1)/2);
y=y-(Geometry.nVoxel(2)/2);
mask=zeros(Geometry.nVoxel(1),Geometry.nVoxel(2));
mask(x.^2+y.^2<=(Geometry.nVoxel(1)/2).^2)=1;
imgFDKpost=bsxfun(@times,mask,resFDK);
imgFDKpost(imgFDKpost<0)=0;
% imgFDKpost=smooth3(imgFDKpost,'gaussian');
%%

% figure(2);imshow(imgFDKoff(:,:,ii),[0 0.05]);
%  plotImg(imgFDK,1,'Z');

%% SART
% 
% 
% tic
% [res,err]=SART_CBCT(data,Geometry,alpha,60);
% toc