function plotImg(img,step,cross)

if nargin<3
    cross=1;
end

if cross==3 || strcmp(cross,'Z')
    img=permute(img,[3 2 1]);
end

for ii=size(img,1):-1*step:1
    imagesc((squeeze(img(ii,:,:)))'); 
    
    
    axis image; 
    axis equal; 
    
    colormap(py_B_cmap); 
    colorbar; 
    caxis([min(img(:)),max(img(:))]);
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
end