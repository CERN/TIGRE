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
Geometry.accuracy=1;





 alpha=rand(1,20)*2*pi;
 
 % warm up
 img=ones(Geometry.nVoxel');
 b=Ax(img,Geometry,alpha);
 dump=Atb(b,Geometry,alpha);
 for ii=0:4
     Geometry.nDetector=[64; 64]*2^ii;
     Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;
     for jj=0:4

         Geometry.nVoxel=[64;64;64]*2^jj;
         Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;
         img=rand(Geometry.nVoxel');       
         
         b=Ax(img,Geometry,alpha);
         tic;
         b=Ax(img,Geometry,alpha);
         t(ii+1,jj+1)=toc/length(alpha);
         
         clear img;
         
         dump=Atb(b,Geometry,alpha);
         tic;
         dump=Atb(b,Geometry,alpha);
         tb(ii+1,jj+1)=toc/length(alpha);
         clear dump b;
     end
 end
 disp('phase 1 ended');
 
 Geometry.accuracy=0.1;
  for ii=0:4
     Geometry.nDetector=[64; 64]*2^ii;
     Geometry.dDetector=Geometry.sDetector./Geometry.nDetector;
     for jj=0:4
         Geometry.nVoxel=[64;64;64]*2^jj;
         Geometry.dVoxel=Geometry.sVoxel./Geometry.nVoxel;
         img=rand(Geometry.nVoxel');
         
         b=Ax(img,Geometry,alpha);
         tic;
         b=Ax(img,Geometry,alpha);
         t2(ii+1,jj+1)=toc/length(alpha);
         
         
         clear img;

         dump=Atb(b,Geometry,alpha);
         tic;
         dump=Atb(b,Geometry,alpha);
         tb2(ii+1,jj+1)=toc/length(alpha);
         clear dump b;
     end
  end
 %% conver to ms
  t=t*1000;
  t2=t2*1000;
  tb=tb*1000;
  tb2=tb2*1000;
  break
 %%
 tplot1=tb';
 tplot2=tb2';
 addpath C:\Users\aabhca20\Documents\perceptually_uniform_colormaps

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])

 subplot(121)
 
% bar3(t);
% zlim([0.74 10000])
% set(gca,'Zscale','log')
% pause(1)
% ticks=get(gca,'Ztick');
% ticklabel=(get(gca,'ZtickLabel'));
% set(gca,'Zscale','linear')
% cla
%  
 cmap=py_C_cmap();
%  cmap=viridis();
%  cmap='jet';
cmap=[cmap; repmat(cmap(end,:),[100,1])];
 
 b=bar3(log10(tplot1));
 zlim([0 log10(15000)])
 set(gca,'Ztick',log10([1,10,100,1000,10000, 12000]));
 set(gca,'ZtickLabel',{'1ms','10ms','100ms','1s','10s',''});
 title('Low accuracy')
 xlabel('Detector size')
 grid on
 ax=gca;
 ax.XTickLabel={'64^2','128^2','256^2','512^2','1024^2'};
 ax.YTickLabel={'64^3','128^3','256^3','512^3','1024^3'};
 ylabel('Image size')
 zlabel('Time')
 colormap(cmap);
 ax.CLim=[0 log10(10000)];
 for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
 end
 view(-133,26)

 subplot(122)
 

 
 
 b=bar3(log10(tplot2));
 zlim([0 log10(12000)])
 set(gca,'Ztick',log10([1,10,100,1000,10000,12000]));
 set(gca,'ZtickLabel',{'1ms','10ms','100ms','1s','','12s'});
 title('High accuracy')
 xlabel('Detector size')
 grid on
 ax=gca;
 ax.XTickLabel={'64^2','128^2','256^2','512^2','1024^2'};
 ax.YTickLabel={'64^3','128^3','256^3','512^3','1024^3'};
 ylabel('Image size')
 zlabel('Time')
 colormap(cmap);
 ax.CLim=[0 log10(10000)];
 for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
 end
 view(-133,26)