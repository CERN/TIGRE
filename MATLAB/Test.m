


%% Apply reconstruction full circle Traj Alderson %%%% 


load FullSinglescanAlderson 
load FullSinglescanAldersonangle 
load FullSinglescanAldersonangleRad

 
%%% apply inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxdata=max(FullSinglescanAlderson (:)); % =1023
    
for ii=1:size(FullSinglescanAlderson,3)

FullSinglescanAldersoninverse(:,:,ii)=(-1*single(FullSinglescanAlderson(:,:,ii)))+255;
end

figure,imshow(FullSinglescanAldersoninverse(:,:,36),[])




geo.DSD =1195;                             % Distance Source Detector      (mm)
%geo.DSD =895;
geo.DSO = 810;                             % Distance Source Origin        (mm)
geo.nDetector=[512; 384];
%geo.nDetector=[512; 512];

% number of pixels              (px)
%geo.dDetector=[0.7413; 0.7413];
geo.dDetector=[0.7420; 0.7420];

geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[256;256;256]; %the best  results               % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)


geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
geo.offOrigin=[0;0;0];
  





  clear projections anglesall
  
%%% Best 1 %%%%  
projections=FullSinglescanAldersoninverse;
anglesall=FullSinglescanAldersonangleRad;
% plotProj([ projections(:,:,40)],anglesall); 
 projections([1:64 449:512],:,:)=[];  %%%%for cutting detector parts
% projections([1:64 449:512],:,:)=58;  %%%% better 
% 
projections=projections(:,:,1:5:end);
anglesall=anglesall(:,1:5:end);


imgOSSART=OS_SART(projections,geo,anglesall,60);  % 30 iter, 3:35 % 60 iter 6:30s
imgSIRT=SIRT(projections,geo,anglesall,30);       % 30 iter, 0:32
imgSIRT150=SIRT(projections,geo,anglesall,150);   %150 iter, 3:48
imfFDK=FDK(projections,geo,anglesall);            % "free"
[imgLSQR,res1]=LSQR(projections,geo,anglesall,60);% 30 iter, 0:30
[imgCGLS,res2]=CGLS(projections,geo,anglesall,30);% 30 iter  0:49

[imgLSMR,res3]=LSMR(projections,geo,anglesall,30);% 30 iter  0:30
[imgLSMR2,res4]=LSMR(projections,geo,anglesall,30,'lambda',30);% 30 iter  0:30
[imghLSQR, res5]=hybrid_LSQR(projections,geo,anglesall,30, 'lambda', 10); % 30 iter  0:40

%% clean up
imgOSSART=cropCBCT(imgOSSART(:,:,30:end-29),geo);
imgSIRT=cropCBCT(imgSIRT(:,:,30:end-29),geo);
imgSIRT150=cropCBCT(imgSIRT150(:,:,30:end-29),geo);
imfFDK=cropCBCT(imfFDK(:,:,30:end-29),geo);
imgLSQR=cropCBCT(imgLSQR(:,:,30:end-29),geo);
imgCGLS=cropCBCT(imgCGLS(:,:,30:end-29),geo);
imgLSMR=cropCBCT(imgLSMR(:,:,30:end-29),geo);
imgLSMR2=cropCBCT(imgLSMR2(:,:,30:end-29),geo);
imghLSQR=cropCBCT(imghLSQR(:,:,30:end-29),geo);
% imfFDK=imfFDK(:,:,30:end-29);

%
zer=zeros(size(imfFDK),'single');
%%
plotImg([imgLSMR,imgLSMR2,imghLSQR;imgOSSART,imgCGLS, imgLSQR;imfFDK, imgSIRT,imgSIRT150;],'dim',3,'clims',[0 3],'slice',205-30)  
plotImg(cat(3,flip([imgLSMR; imgLSMR2; imghLSQR],3),...
              flip([imgOSSART;imgCGLS;imgLSQR],3),...
              flip([imfFDK;imgSIRT;imgSIRT150],3)...
              ),'dim',2,'clims',[0 3],'slice',128);
%%
return
 %% TV experiments
 
[imgIRN_TV_CGLS,res1]=IRN_TV_CGLS(projections,geo,anglesall,60,'lambda',50,'niter_outer',4);
imgIRN_TV_CGLS=cropCBCT(imgIRN_TV_CGLS(:,:,30:end-29),geo);
[imgIRN_TV_CGLS2,res2]=IRN_TV_CGLS(projections,geo,anglesall,60,'lambda',20,'niter_outer',4);
imgIRN_TV_CGLS2=cropCBCT(imgIRN_TV_CGLS2(:,:,30:end-29),geo);
[imgIRN_TV_CGLS3,res3]=IRN_TV_CGLS(projections,geo,anglesall,60,'lambda',5,'niter_outer',4);
imgIRN_TV_CGLS3=cropCBCT(imgIRN_TV_CGLS3(:,:,30:end-29),geo);
imgOSASDPOCS=OS_ASD_POCS(projections,geo,anglesall,60); %OSC options
imgOSASDPOCS=cropCBCT(imgOSASDPOCS(:,:,30:end-29),geo);

% imghLSQR=hybrid_fLSQR_TV(projections,geo,anglesall,30,'lambda',50);
% imghLSQR=cropCBCT(imghLSQR(:,:,30:end-29),geo);
