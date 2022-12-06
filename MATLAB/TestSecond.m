
load fullcircallspeed 
load fullcircallanglespeed 
load fullcircallangleRadspeed 



maxdata=max(fullcircallspeed (:)); % =1023
    
for ii=1:size(fullcircallspeed ,3)

fullcircallspeedinverse(:,:,ii)=(-1*single(fullcircallspeed (:,:,ii)))+255;
end

figure,imshow(fullcircallspeedinverse(:,:,136),[])

    
clear projections anglesall
  
%%% Best 1  158 angles %%%%  

projections=fullcircallspeedinverse;
anglesall=fullcircallangleRadspeed;

%projections=fullcircallspeedinverse(:,:,[1:2:329 1:4:300]);
%anglesall=fullcircallangleRadspeed(:,[1:2:329 1:4:300]);
% plotProj([ projections(:,:,40)],anglesall); 
 projections([1:64 449:512],:,:)=[];  %%%%for cutting detector parts
 %projections([1:64 449:512],:,:)=43;  %%%% better 

 

geo.DSD =1195;                             % Distance Source Detector      (mm)
%geo.DSD =895;
geo.DSO = 810;                             % Distance Source Origin        (mm)
geo.nDetector=[512; 384];
%geo.nDetector=[512; 512];

% number of pixels              (px)
geo.dDetector=[0.7413; 0.7413];
%geo.dDetector=[0.7420; 0.7420];

geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[256;256;256]; %the best  results               % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)


geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
geo.offOrigin=[0;0;0];
  
 
lambda=1;
lambdared=0.98;
verb=true;verbose=true
TVlamb=600; %%% ok value 

epsilon=10^6;
alpha=0.002;
ng=25;
ratio=0.94;

 ASDPOCSRealdataTHoraxcircular240range=ASD_POCS(single(projections),geo,anglesall,10,...
    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb); % less important.
