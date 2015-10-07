function plotProj(proj,alpha)

for ii=1:size(proj,3)
    image=squeeze(proj(:,:,ii));
    figure(1); imagesc((image));axis image;  axis equal; colormap gray; colorbar; xlabel('-> U');ylabel('-> V');set(gca,'XTick',[]);set(gca,'YTick',[]);set(gca,'YDir','normal');
    title(['Degree : ',num2str(alpha(ii)*180/pi)]);
    pause(0.01);
end
