function plotProj(proj,alpha)
figure();
for ii=1:size(proj,3)
    image=squeeze(proj(:,:,ii));
     imagesc((image));
    
    axis image;  
    axis equal; 
    
    colormap(py_B_cmap); 
    colorbar; 
    
    xlabel('-> U');
    ylabel('-> V');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YDir','normal');
    
    title(['Degree : ',num2str(alpha(ii)*180/pi)]);
    pause(0.01);
end
