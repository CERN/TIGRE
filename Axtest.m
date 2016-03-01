%% How to use Ax with CUDA:
%
%
% Attention: The code has been written for a NVIDIA Tesla k40 graphic card,
% and most likely need to be tunned if it is wanted to use with another a
% different one.
%
%
% This code solves Geometrymetrically the A*x of a linear system of equations
% for the Cone Beam Computed Tomography Geometrymetry. In other words, it
% computes a bection with a given image. It can be used for both forward
% bection (simulating the measured 2D bections) and for iterative
% inverse solving of the image, where in many cases (as in SART algorithm)
% the A*x term has to be computed several times. The benefit is that there
% is no need of actually computing the matrix A.
%
% The code makes use of graphic cards hardware implemented texture memory
% and the multy CPU arquitechture of modern GPUs to simultaneously compute
% the bection in each detector pixel.
%
%
% A single call to Ax will compute a single bection (as for 05/08/2015).
% Detector array is plane.
%
% INPUT ARGUMETS:
%       bection=Ax(image, Geometrymetry)
%
% IMAGE: a 3D image of the current solution.
% GeometryMETRY: an structure describing the Geometrymetry of the CBCT, containing:
%           -nVoxel:        3x1 matrix of number of voxels in the image
%           -sVoxel:        3x1 matrix with the total size in mm of the image
%           -dVoxel:        3x1 matrix with the size of each of the voxels in mm
%           -nDetector:     2x1 matrix of number of voxels in the detector plane
%           -sDetector:     2x1 matrix with the total size in mm of the detector
%           -dDetector:     2x1 matrix with the size of each of the pixels in the detector in mm
%           -DSD:           1x1 scalar value. Distance Source Detector, in mm
%           -DSO:           1x1 scalar value. Distance Source Origin.
%           -offOrigin:     3x1 matrix with the offset in mm of the centre of the image from the origin.
%           -offDetector:   2x1 matrix with the offset in mm of the centre of the detector from the x axis
% ALPHA: the angle in radians of the bection
%

%%  GeometryMETRY DEFINITION
%
%                  Detector plane, behind
%              |-----------------------------|
%              |                             |
%              |                             |
%              |                             |
%  Centered    |                             |
%    at O      A V    +--------+             |
%              |     /        /|             |
%     A Z      |    /        / |*D           |
%     |        |   +--------+  |             |
%     |        |   |        |  |             |
%     |        |   |     *O |  +             |
%     *--->y   |   |        | /              |
%    /         |   |        |/               |
%   V X        |   +--------+        U       |
%              .--------------------->-------|
%
%            *S
%
%
%
%
%
%
% %% Example of use
close all
clear
InitToolbox;
%%
% testtime=0;


%
% Geometrymetry.nVoxel=[128;128;128];
% Geometrymetry.sVoxel=[460;460;460];
% Geometrymetry.dVoxel=Geometrymetry.sVoxel./Geometrymetry.nVoxel;
%
% Geometrymetry.nDetector=[128;128];
% Geometrymetry.sDetector=[ 1024;800];
% Geometrymetry.dDetector=Geometrymetry.sDetector./Geometrymetry.nDetector;
%
% Geometrymetry.DSD = 1500;
% Geometrymetry.DSO = 1100;
%
% Geometrymetry.offOrigin=[0; 0; 0];
% Geometrymetry.offDetector=[0;0];
% Geometrymetry.accuracy=1;

%% P from matrix code?
% clear Geometrymetry

%% Geometry structure definition.
% Distances
geo.DSD = 1536;
geo.DSO = 1000;
% Detector parameters
geo.nDetector=[512; 512];					% number of pixels 
geo.dDetector=[0.8; 0.8]; 					% size in mm of each pixel
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector in mm
% Image parameters
geo.nVoxel=[256;256;256];                   % number of voxels in the image
geo.sVoxel=[256;256;256];                   % total size of the image in mm
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size in mm of each voxel
% Offsets
geo.offOrigin =[0;0;0];                     % V_orig
geo.offDetector=[0; 0];                     % V_det
Geometry.accuracy=0.1;
%
% [P,~] = xread('C:\VOL_CT_modified\rando_head\');
% alpha=
% % %% Real image in the coords we like
load img128
img=double(img);

[y, x, z]=...
    ndgrid(linspace(1,size(img,1),Geometry.nVoxel(1)),...
    linspace(1,size(img,2),Geometry.nVoxel(2)),...
    linspace(1,size(img,3),Geometry.nVoxel(3)));
imOut=interp3(img,x,y,z);
img=imOut;
% img=phantom3dAniso(Geometry.nVoxel);
clear imgOut x y z
%% plot image
%  plotImg(img,5)

% img=ones(Geometrymetry.nVoxel');

%  img=ones(Geometrymetry.nVoxel');
%  alpha=-pi/2;

%% bect

alpha=[0:20:359]*pi/180;
% % img=phantom3dAniso(Geometrymetry.nVoxel);
% % img=img./max(img(:));
% % tic
b=Ax(img,Geometry,alpha,'Krylov');
% % toc
% % b=b+(randn(size(b))-0.5)*max(b(:))/80;
b=addCTnoise(b);
% % 
% 
% 
% 
% 
% sizeZ=Geometry.sVoxel(3)*2/3+Geometry.sVoxel(3)/5*sin(alpha);

% tic
[imgOSSART]=OS_SART_CBCT(b,Geometry,alpha,40);
[imgfdk]=FDK_CBCT(b,Geometry,alpha);
[imgADSPOCS]=ADS_POCS_CBCT(b,Geometry,alpha,50,im3Dnorm(imgOSSART,'L2'));
[imgBADSPOCS]=B_ADS_POCS_beta_CBCT(b,Geometry,alpha,10,5,im3Dnorm(imgOSSART,'L2'),0.75);

break;
[imgCGLS,errCGLS]=CGLS_CBCT(b,Geometrymetry,alpha,9);

% [imgADSPOCS,tv]=ADS_POCS_CBCT(b,Geometrymetry,alpha,50,im3Dnorm(imgOSSART,'L2'));
% [imgBADSPOCS,tv]=B_ADS_POCS_beta_CBCT(b,Geometrymetry,alpha,50,im3Dnorm(imgOSSART,'L2'),0.75);

% imshow([imgADSPOCS(:,:,60) imgOSSART(:,:,60)],[])
%%
% %
% s=512;
%
% f0=ones(s,s,s)*1;%*rand(1);%/s^3;
% f0(5:100,5:100,5:100)=0;
% % f0=imgOSSART;
% dtvg=6;
% ng=100;
% f=f0;
% % tic
% % for ii=1:ng
% %     % Steepest descend of TV norm
% %     df=gradientTVnorm(f,'backward');
% % %     n(ii)=im3Dnorm(df,'L2')
% %     df=df./im3Dnorm(df,'L2');
% %     f=f-dtvg.*df;
% % end
% % toc
% tic
% pocsed=minimizeTV(f0,dtvg,ng);
% toc
%
% % close all
% % imagesc(f(:,:,10)); colormap gray; axis xy
% % figure
% imagesc(pocsed(:,:,10)); colormap gray;axis xy

break
%%
tic
[imgCGLS,errCGLS]=CGLS_CBCT(b,Geometrymetry,alpha,8);
toc
tic
[imgSART,errSART]=OS_SART_CBCT(b,Geometrymetry,alpha,30,'BlockSize',20,'Init','multigrid');
toc
% break;
break;

tic
[imgSART,errSART]=OS_SART_CBCT(b1,Geometrymetry,alpha,30,'BlockSize',20);
toc
%  tic;
%  b2=Ax(img,Geometrymetry,alpha);
%  toc;
break

break
maxb=max(b(:));
% bnoise=imnoise(b/maxb,'poisson');
% bnoise=bnoise.*maxb;
% b=bnoise;
% %
%
%  %FDK
%  tic
% Geometrymetry.filter='ram-lak';
% b_filt = filtering(b,Geometrymetry,alpha); % Not sure if offsets are good in here
% Geometrymetry=rmfield(Geometrymetry,'filter');
% imgFDK=Atb(b_filt,Geometrymetry,alpha);
% toc

tic
[imgCGLS,errCGLS]=CGLS_CBCT(b,Geometrymetry,alpha,60);
toc
% tic
% [imgSART,errSART]=SART_CBCT(b,Geometrymetry,alpha,30);
% toc

break
plotb(b,alpha);



%% ALGORITMS!
%
tic
[res,err]=CGLS_CBCT(b,Geometrymetry,alpha,60);
toc
tic
[res,err]=SART_CBCT(b,Geometrymetry,alpha,10);
toc

plot(err);
break
plotImg(res,1,3);
plotImg(img,1,3);
