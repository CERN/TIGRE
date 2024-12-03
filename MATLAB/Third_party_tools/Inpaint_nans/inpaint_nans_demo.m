%% Surface Fit Artifact Removal

%% Construct the Surface
[x,y] = meshgrid(0:.01:1);
z0 = exp(x+y);

close all
figure
surf(z0)
title 'Original surface'

znan = z0;
znan(20:50,40:70) = NaN;
znan(30:90,5:10) = NaN;
znan(70:75,40:90) = NaN;

figure
surf(znan)
title 'Artifacts (large holes) in surface'

%% In-paint Over NaNs
z = inpaint_nans(znan,3);
figure
surf(z)
title 'Inpainted surface'

figure
surf(z-z0)
title 'Inpainting error surface (Note z-axis scale)'

%% Comapre to GRIDDATA
k = isnan(znan);
zk = griddata(x(~k),y(~k),z(~k),x(k),y(k));
zg = znan;
zg(k) = zk;

figure
surf(zg)
title(['Griddata inpainting (',num2str(sum(isnan(zg(:)))),' NaNs remain)'])

figure
surf(zg-z0)
title 'Griddata error surface'
