%% Demo 23: 4D reconstruction
%
% 4D CT happens with both introduction of motion, or spectral imaging. 
% This demo shows how to simulate and recosntruct 4D CT using TIGRE.
% Currently TIGRE only contains a single 4D algorithm, more to come. Feel
% free to add more!
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------
%% Initialization
clear all

geo=defaultGeometry('nVoxel',[128,128,128]','nDetector',[128,128]);
angles=linspace(0,2*pi,50);

%% lets simulate some 4D data.
% While it makes no anatomical sense, lets make the bones in the body to
% change slightly with e.g. X-ray spectra 
% (all tissues are energy dependent, we are just making up data and values
% here)
phan=headPhantom(geo.nVoxel);
mask=phan>0.5;
phantom=repmat(phan,[1,1,1,5]);
phantom(:,:,:,2)=phan.*~mask+phan.*mask.*1.05;
phantom(:,:,:,3)=phan.*~mask+phan.*mask.*1.1;
phantom(:,:,:,4)=phan.*~mask+phan.*mask.*1.2;
phantom(:,:,:,5)=phan.*~mask+phan.*mask.*1.25;

% Lets forward project. 
for ii=1:size(phantom,4)
    data(:,:,:,ii)=Ax(phantom(:,:,:,ii),geo,angles);
end
N       = [geo.nVoxel' size(data,4)]; 



% Recover the image
mu      = 1;
lambda  = 1;
gamma   = .1;
alpha   = [1 1 1];  % aloha(1) for x and y, alpha(2) for z and alpha(3) for time (or freuqency)
nInner  = 2;
nBreg   = 5;

[split_bregg] = extended_split_bregman_TV(data,geo,angles,mu,lambda,gamma,alpha,nInner,nBreg,1e-4,100);

%% Recosntruct with FDK
for ii=1:size(phantom,4)
    imgFDK(:,:,:,ii)=FDK(data(:,:,:,ii),geo,angles);
end
for ii=1:size(phantom,4)
    imgOSSART_TV(:,:,:,ii)=OS_ASD_POCS(data(:,:,:,ii),geo,angles,20,'maxL2err',0.2*im3Dnorm(Ax(imgFDK(:,:,:,ii),geo,angles)-data(:,:,:,ii),'L2'));
end
%% Compute error images
error_TV=abs(split_bregg-phantom);
error_FDK=abs(imgFDK-phantom);
error_OSSART_TV=abs(imgOSSART_TV-phantom;

for ii=1:size(phantom,4)
    UQI_TV(ii)=UQI(split_bregg(:,:,:,ii),phantom(:,:,:,ii));
    UQI_FDK(ii)=UQI(imgFDK(:,:,:,ii),phantom(:,:,:,ii));
    UQI_OSSART_TV(ii)=UQI(imgOSSART_TV(:,:,:,ii),phantom(:,:,:,ii));
end
%% Plot
plotImg( [split_bregg(:,:,:,1),imgFDK(:,:,:,1), imgOSSART_TV(:,:,:,1); ...
          split_bregg(:,:,:,3),imgFDK(:,:,:,3), imgOSSART_TV(:,:,:,3); ...
          split_bregg(:,:,:,5),imgFDK(:,:,:,5), imgOSSART_TV(:,:,:,5) ],'Dim','z','slice',34);
plotImg( [error_TV(:,:,:,1),error_FDK(:,:,:,1),error_OSSART_TV(:,:,:,1);...
          error_TV(:,:,:,3),error_FDK(:,:,:,3),error_OSSART_TV(:,:,:,3);...
          error_TV(:,:,:,5),error_FDK(:,:,:,5),error_OSSART_TV(:,:,:,5) ],'Dim','z','slice',34);
