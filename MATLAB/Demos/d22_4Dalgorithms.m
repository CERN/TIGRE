clear all

geo=defaultGeometry('nVoxel',.25*[512,512,512]','nDetector',.25*[512,512]);

angles=linspace(0,2*pi,50);
phan=headPhantom(geo.nVoxel);
%phan=thoraxPhantom(geo.nVoxel);
proj=Ax(phan,geo,angles,'interpolated');

% return
data(:,:,:,1)=proj;
data(:,:,:,2)=proj;
data(:,:,:,3)=proj;
data(:,:,:,4)=proj;
data(:,:,:,5)=proj;
N       = [geo.nVoxel' size(data,4)]; 



% Recover the image
mu      = 1;
lambda  = 1;
gamma   = .1;
alpha   = [1 1 1];  % aloha(1) for x and y, alpha(2) for z and alpha(3) for time (or freuqency)
nInner  = 1;
nBreg   = 3;

[u] = extended_split_bregman_TV(data,geo,angles, N,mu,lambda,gamma,alpha,nInner,nBreg);

imgFDK=FDK(data(:,:,:,1),geo,angles);


for ii=1:size(data,4)
    uu=u(:,:,:,ii);
%     slice(uu,[fix(geo.nVoxel/2)],[fix(geo.nVoxel/2)],fix(geo.nVoxel/2));shading interp;colormap jet;pause(0.50)
%     drawnow
end
plotImg( [u(:,:,:,1);imgFDK ],'Dim','z','slice',64);
