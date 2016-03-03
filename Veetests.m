clear; 
clc;

tic
InitToolbox;

Geometry.nVoxel=[128;128;128];                              
Geometry.sVoxel=[460;460;460];                             
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;


Geometry.nDetector=[128;128];                               
Geometry.sDetector=[ 1024;800];                             
Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;  
Geometry.DSD = 1500;  
Geometry.DSO = 1100; 

Geometry.offOrigin=[0; 0; 0];           
Geometry.offDetector=[0;0];       
Geometry.accuracy=1;  


alpha=[0:1:359]*pi/180;



phantom=thoraxPhantom(Geometry.nVoxel);
% phantom=ones(Geometry.nVoxel');
% phantom(80:100,80:100,80:100)=0;
% phantom2=phantom*0.99;

%Generate data
data=Ax(phantom,Geometry,alpha); % this uses GPU




% plot the projections jumping 1 and save the result as a gif

% plotProj(data,alpha)
%% Reconstruct image

% FDK
%Input is data,geosize,angles
% [resFDK,errFDK,rmseFDK,ccFDK,mssimFDK]=FDK_CBCT(data,Geometry,alpha,phantom);
% 
%  u=UQI(phantom,resFDK);
% 


% niter=200;
niter=20;
nblock=20;

% [resSIRT,errorSIRT,rmseSIRT,ccSIRT,mssimSIRT]=SIRT_CBCT(data,Geometry,alpha,niter);
[resOSSART,errorOSSART,rmseOSSART,ccOSSART,mssimOSSART,uqiOSSART]=OS_SART_CBCT(data,Geometry,alpha,niter,'BlockSize',nblock);
% [resSART,errorSART,rmseSART]=SART_CBCT(data,Geometry,alpha,niter);

toc
break
%Plot the RMSE
figure
plot(1:niter,rmseOSSART);xlabel('Number of Iteration');ylabel('RMSE');title('OS-SART with 180 projections')

%Plot the Pearson correlation coefficient (CC)
figure
plot(1:niter,ccOSSART);xlabel('Number of Iteration');ylabel('Correlation coefficient');%title('OS-SART')

%Plot the mean structural similarity index (MSSIM)
figure
plot(1:niter,mssimOSSART);xlabel('Number of Iteration');ylabel('Mean structural similarity index')

break


rmFDK=RMSE(phantom,resFDK,Geometry);
rmOSSART=RMSE(phantom,resOSSART,Geometry);


figure
plot(1:Geometry.nVoxel(1)*Geometry.nVoxel(2)*Geometry.nVoxel(3),rmFDK,1:Geometry.nVoxel(1)*Geometry.nVoxel(2)*Geometry.nVoxel(3),rmOSSART)

%Plot RMSE VS Cycles
%Image row plots: The plot running through an image at specific x,y or z
%axis


break






%% Plot the result image

% plot XY slices
plotImg(resFDK,'Dim','Z')

% plotImg(resSIRT,'Dim','X')