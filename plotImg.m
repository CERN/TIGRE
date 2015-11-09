function plotImg(img,step,cross,savegif)

if nargin<3
    cross=1;
    savegif=0;
end
if nargin <4
    savegif=0;
end
if cross==3 || strcmp(cross,'Z')
    img=permute(img,[3 2 1]);
end
climits=prctile(img(:),[1 99]);

fh=figure();
filename='sart30.gif';
for ii=size(img,1):-1*step:1
    imagesc((squeeze(img(ii,:,:)))'); 
    
    
    axis image; 
    axis equal; 
    
    colormap('gray'); 
    colorbar; 
    caxis([climits(1),climits(2)]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YDir','normal');
    
    
    if cross==3 || strcmp(cross,'Z')
        xlabel('->Y');
        ylabel('X<-');
        title(['Top to bottom ->Z : ',num2str(ii)]);
    else
        xlabel('->Y');
        ylabel('->Z');
        title(['Source to Detector direction ->X : ',num2str(ii)]);
    end
    drawnow
    
    if savegif
        
      frame = getframe(fh);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if ii == 256;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
      end
    end
end

