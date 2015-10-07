function plotImg(img,step)

for ii=size(img,1):-1*step:1
    imagesc((squeeze(img(ii,:,:)))'); axis image; axis equal; colormap gray; colorbar; caxis([min(img(:)),max(img(:))]);set(gca,'XTick',[]);set(gca,'YTick',[]);set(gca,'YDir','normal');
    xlabel('->Y');
    ylabel('->Z');
    title(['Source to Detector direction ->X : ',num2str(ii)]);
    drawnow
end