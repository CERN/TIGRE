%% Tests the performance of the projection adn backprojection codes.
%
%
%
% IMPORTANT NOTE: In order to be able to run this code and get results, you
% need to go to the 3 .cu files that do the projection adn backprojection
% work ,"ray_interpolated_projection.cu" , "Siddon_projection.cu", and
% ,"voxel_backprojection.cu" and in each of
% them, set the variable "timethis" to true. Then, recompile, and run this
% code.
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
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------

Geometry.DSD = 1536;   
Geometry.DSO = 1000;

Geometry.nDetector=[512; 512];
Geometry.dDetector=[0.8; 0.8];
Geometry.sDetector=Geometry.nDetector.*Geometry.dDetector;

Geometry.nVoxel=[512;512;512];
Geometry.sVoxel=[256;256;256];
Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;

Geometry.offOrigin=[0;0;0]; 
Geometry.offDetector=[0;0];
Geometry.accuracy=0.5;





 alpha=linspace(0,pi,100); 
 alpha=alpha(1:32);
%  for ii=0:4
%      Geometry.nDetector=[64; 64]*2^ii;
%      Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;
%      for jj=0:4
%          Geometry.nVoxel=[64;64;64]*2^jj;
%          Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;
%          img=single(rand(Geometry.nVoxel')); 
% %          tic
% %          b=Ax(img,Geometry,alpha);
% %           t_interpolation(ii+1,jj+1)=toc;
%          t_interpolation(ii+1,jj+1)=str2double(evalc('b=Ax(img,Geometry,alpha);'));
%          t_backprojection(ii+1,jj+1)=str2double(evalc('Atb(b,Geometry,alpha);'));
% %          tic
% %          Atb(b,Geometry,alpha);
% %           t_backprojection(ii+1,jj+1)=toc;
%          clear  b;
%      end
%  end
%  
 Geometry.accuracy=0.5;
 for ii=0:4
     Geometry.nDetector=[64; 64]*2^ii;
     Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;
     for jj=0:4
         Geometry.nVoxel=[64;64;64]*2^jj;
         Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;
         img=single(rand(Geometry.nVoxel'));       
         t_siddon(ii+1,jj+1)=str2double(evalc('b=Ax(img,Geometry,alpha,''ray-voxel'');'));
         t_backprojection2(ii+1,jj+1)=str2double(evalc('Atb(b,Geometry,alpha,''matched'');'));
         t_backprojection(ii+1,jj+1)=str2double(evalc('Atb(b,Geometry,alpha,''FDK'');'));

%           t_siddon(ii+1,jj+1)=toc;
         clear  b;
     end
 end
 break
 %% Plot projection and bakcprojection performance.
 tplot1=t_interpolation'*100; % base unit 1ms
 tplot2=t_siddon'*100;
 tplot3=t_backprojection';

 
 
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])


subplot(121)
 

 cmap=viridis();
cmap=[cmap; repmat(cmap(end,:),[100,1])];
 
 b=bar3(log10(tplot1));
 zlim([0 log10(20000)])
 set(gca,'Ztick',log10([1,10,100,1000,10000,20000]));
 set(gca,'ZtickLabel',{'10ns','100ns','1ms','10ms','100ms','200ms'});
 title('Projection with Interpolation','fontsize',25)
 xlabel('Detector size','fontsize',25)
 grid on
 ax=gca;
  set(ax,'fontsize',20);

 ax.XTickLabel={'64^2','128^2','256^2','512^2','1024^2'};
 ax.YTickLabel={'64^3','128^3','256^3','512^3','1024^3'};
 ylabel('Image size','fontsize',25)
 zlabel('Time','fontsize',25)
 colormap(cmap);
 ax.CLim=[0 log10(10000)];
 for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
 end
 set(gca,'DataAspectRatio',[1 1 0.9]);
 xlim([0,5.4])
 ylim([0,5.4])
 view(-133,26)

 subplot(122)
 
 b=bar3(log10(tplot2));
  zlim([0 log10(20000)])
 set(gca,'Ztick',log10([1,10,100,1000,10000,20000]));
 set(gca,'ZtickLabel',{'10ns','100ns','1ms','10ms','100ms','200ms'});
 title('Projection with ray-voxel intersection','fontsize',25)
 xlabel('Detector size','fontsize',25)
 grid on
 ax=gca;
  set(ax,'fontsize',20);

 ax.XTickLabel={'64^2','128^2','256^2','512^2','1024^2'};
 ax.YTickLabel={'64^3','128^3','256^3','512^3','1024^3'};
 ylabel('Image size','fontsize',25)
 zlabel('Time','fontsize',25)
 colormap(cmap);
 ax.CLim=[0 log10(10000)];
 for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
 end
 set(gca,'DataAspectRatio',[1 1 0.9]);
 xlim([0,5.4])
 ylim([0,5.4])
 view(-133,26)

 %%
 % Backprojection plot.
 tplot3=t_backprojection'*100/32;
 tplot4=t_backprojection2'*100/32;
 
 figure(2)
 subplot(121)
set(gcf,'units','normalized','outerposition',[0 0 1 1])


 cmap=plasma();
 cmap=[cmap; repmat(cmap(end,:),[100,1])];
 b=bar3(log10(tplot3));
 zlim([0 log10(20000)])
 set(gca,'Ztick',log10([1,10,100,1000,10000,20000]));
 set(gca,'ZtickLabel',{'10ns','100ns','1ms','10ms','100ms','200ms'});
 title('Backprojection (FDK)','fontsize',25)
 xlabel('Detector size','fontsize',25)
 grid on
 ax=gca;
 ax.XTickLabel={'64^2','128^2','256^2','512^2','1024^2'};
 ax.YTickLabel={'64^3','128^3','256^3','512^3','1024^3'};
 set(ax,'fontsize',20);
 ylabel('Image size','fontsize',25)
 zlabel('Time','fontsize',25)
 colormap(cmap);
 ax.CLim=[0 log10(20000)];
 for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
 end
 set(gca,'DataAspectRatio',[1 1 1]);
 xlim([0,5.4])
 ylim([0,5.4])
 view(-133,26)
 
  subplot(122)
set(gcf,'units','normalized','outerposition',[0 0 1 1])


 cmap=plasma();
 cmap=[cmap; repmat(cmap(end,:),[100,1])];
 b=bar3(log10(tplot4));
 zlim([0 log10(20000)])
 set(gca,'Ztick',log10([1,10,100,1000,10000,20000]));
 set(gca,'ZtickLabel',{'10ns','100ns','1ms','10ms','100ms','200ms'});
 title('Backprojection (pseudo-matched)','fontsize',25)
 xlabel('Detector size','fontsize',25)
 grid on
 ax=gca;
 ax.XTickLabel={'64^2','128^2','256^2','512^2','1024^2'};
 ax.YTickLabel={'64^3','128^3','256^3','512^3','1024^3'};
 set(ax,'fontsize',20);
 ylabel('Image size','fontsize',25)
 zlabel('Time','fontsize',25)
 colormap(cmap);
 ax.CLim=[0 log10(20000)];
 for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
 end
 set(gca,'DataAspectRatio',[1 1 1]);
 xlim([0,5.4])
 ylim([0,5.4])
 view(-133,26)