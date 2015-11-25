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
Geometry.accuracy=0.5;

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
clear D I0 RealProj Proj


%% visualize projections.
vis=0;
if vis
    plotProj(data,alpha);
end
close all


%% reconstruct OS-SART

niter=5;
Geometry.nVoxel=[128;128;128];
Geometry.sVoxel=[256;256;256];
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;


[resSART,errSART]=SART_CBCT(data,Geometry,alpha,niter);
[resOSSART,errOSSART]=OS_SART_CBCT(data,Geometry,alpha,niter,20);
[resSIRT,errSIRT]=SIRT_CBCT(data,Geometry,alpha,niter);

%% Random crap that I migth want to keep
%%
%%
%%
%% Post process
% [x, y]=meshgrid(1:Geometry.nVoxel(1),1:Geometry.nVoxel(2));
% x=x-(Geometry.nVoxel(1)/2);
% y=y-(Geometry.nVoxel(2)/2);
% mask=zeros(Geometry.nVoxel(1),Geometry.nVoxel(2));
% mask(x.^2+y.^2<=(Geometry.nVoxel(1)/2).^2)=1;
% imgFDKpost=bsxfun(@times,mask,imgFDK);
% imgFDKpost(imgFDKpost<0)=0;
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