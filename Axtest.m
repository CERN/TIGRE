%% How to use Ax with CUDA:
% 
%
% Attention: The code has been written for a NVIDIA Tesla k40 graphic card,
% and most likely need to be tunned if it is wanted to use with another a
% different one. 
%
%
% This code solves geometrically the A*x of a linear system of equations
% for the Cone Beam Computed Tomography geometry. In other words, it
% computes a projection with a given image. It can be used for both forward
% projection (simulating the measured 2D projections) and for iterative
% inverse solving of the image, where in many cases (as in SART algorithm)
% the A*x term has to be computed several times. The benefit is that there
% is no need of actually computing the matrix A.
%
% The code makes use of graphic cards hardware implemented texture memory 
% and the multy CPU arquitechture of modern GPUs to simultaneously compute
% the projection in each detector pixel.
%
%
% A single call to Ax will compute a single projection (as for 05/08/2015).
% Detector array is plane.
%
% INPUT ARGUMETS: 
%       projection=Ax(image, geometry)
%  
% IMAGE: a 3D image of the current solution.
% GEOMETRY: an structure describing the geometry of the CBCT, containing:
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
% ALPHA: the angle in radians of the projection
% 

%%  GEOMETRY DEFINITION
%                 
%                  Detector plane, behind
%              |-----------------------------| 
%              |                             | 
%              |                             | 
%              |                             | 
%  Centered    |                             | 
%    at O      |      +--------+             |
%              |     /        /|             |
%     A Z      |    /        / |*D           |
%     |        |   +--------+  |             |
%     |        |   |        |  |             |
%     |        |   |     *O |  +             |
%     *--->y   |   |        | /              |
%    /         |   |        |/               |
%   V X        |   +--------+                |
%              |-----------------------------|
%     
%            *S
%  
%  
%  
%  
%  
%  
%% Example of use
close all
testtime=0;



Geometry.nVoxel=[128;128;128];
Geometry.sVoxel=[460;460;460]; 

Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;

Geometry.nDetector=[256; 200];
Geometry.sDetector=[1024; 800];
Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;

Geometry.DSD = 1500;   
Geometry.DSO = 1100;

Geometry.offOrigin=[0; 0; 0];
Geometry.offDetector=[0; 0];
Geometry.accuracy=0.5;
%%
load img128
% img=double(img);
% ParamSetting;
img1=ones(Geometry.nVoxel')*5;
% img1(10:20,10:20,10:20)=10;
img1(10:20,80:90,10:20)=0;
% img1(1:128,1,1)=1:128;
% img1(:,128,:)=img(:,64,:);
img=img1;


alpha=[0:360]*pi/180;
% alpha=[0,90]*pi/180+pi/2;
% alpha=pi/4;
% alpha=[0 0]*pi/180;
% alpha=30 *pi/180;
tic
b=Ax(img,Geometry,alpha);
toc

for i=1:numel(alpha)
    image=reshape(b(:,i),Geometry.nDetector(2),Geometry.nDetector(1));
    figure(1); imagesc((image)); axis image; axis equal; colormap gray; colorbar; xlabel('-> U');ylabel('-> V');set(gca,'XTick',[]);set(gca,'YTick',[]);set(gca,'YDir','normal');
    title(['Degree : ',num2str(alpha(i)*180/pi)]);
    pause(0.01);
end
break
%%

btest=zeros(Geometry.nDetector(1),Geometry.nDetector(2),size(b,2));
btest=reshape(b,Geometry.nDetector(2),Geometry.nDetector(1),size(b,2));
% alpha=0;
tic
% btest=ones(1,Geometry.nDetector(1),Geometry.nDetector(2));
% btest(1,10:20,10:20)=3;
x=Atb(btest,Geometry,alpha);
toc
% break;
%%
image=reshape(x,Geometry.nVoxel(1),Geometry.nVoxel(2),Geometry.nVoxel(3));
for ii=Geometry.nVoxel(1):-1:1
    imagesc((squeeze(image(ii,:,:)))); axis image; axis equal; colormap gray; colorbar; caxis([min(x),max(x)]);set(gca,'XTick',[]);set(gca,'YTick',[]);set(gca,'YDir','normal');
    xlabel('->Y');
    ylabel('->Z');
    title(['Source to Detector direction ->X : ',num2str(ii)]);
    drawnow
end
% break

%% Test to check if
% b1=A*x1;
% x2=A.'*b1;
% b2=A*x2;
% sum(b1-b2)==0; !!!
% 
close all

alpha=[0 0]*pi/180;
x1=ones(Geometry.nVoxel')*5;
x1(1:30,1:10,1:10)=0;

b1=Ax(x1,Geometry,alpha);

b1=reshape(b1,Geometry.nDetector(2),Geometry.nDetector(1),size(b1,2));
Geometry.offOrigin=[0;0;20];
x2=Atb(b1,Geometry,alpha);

% I need to reshape and transpose x2 now. All this should be avoided, I
% need to change the data input/output system to match
x2=reshape(x2,Geometry.nVoxel(1),Geometry.nVoxel(2),Geometry.nVoxel(3));

b2=Ax(x2,Geometry,alpha);
b2=reshape(b2,Geometry.nDetector(2),Geometry.nDetector(1),size(b2,2));

sum(b1(:)-b2(:));
isequal(b1,b2)

figure
subplot(121)
imagesc(squeeze(b1(:,:,1)));colormap gray;
subplot(122)
imagesc(squeeze(b2(:,:,1)));colormap gray;
break
%% TEst to check backprojection

x=273.8438;
y=-68.4609;
z=-68.4609;

sx=1100;
sy=0;
sz=0;

 vectX=(x -sx); 
 vectY=(y -sy); 
vectZ=(z -sz); 

t=(-400-sx)/vectX;
yres=vectY*t+sy;
zres=vectZ*t+sz;

plot3(x,y,z,'ro')
hold on
plot3(sx,sy,sz,'go')
plot3(-400,yres,zres,'ko')
axis equal
plot3([sx -400],[sy yres],[sz zres],'b')