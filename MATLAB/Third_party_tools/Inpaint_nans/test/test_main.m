%% Repair to an image with 50% random artifacts

% Garden at Sainte-Adresse (Monet, 1867)
garden = imread('monet_adresse.jpg');
G = double(garden);
G(rand(size(G))<0.50) = NaN;
Gnan = G;

G(:,:,1) = inpaint_nans(G(:,:,1),2);
G(:,:,2) = inpaint_nans(G(:,:,2),2);
G(:,:,3) = inpaint_nans(G(:,:,3),2);

figure
subplot(1,3,1)
image(garden)
title 'Garden at Sainte-Adresse (Monet)'

subplot(1,3,2)
image(uint8(Gnan))
title 'Corrupted - 50%'

subplot(1,3,3)
image(uint8(G))
title 'Inpainted Garden'

%% Surface fit artifact removal

[x,y] = meshgrid(0:.01:1);
z0 = exp(x+y);

znan = z0;
znan(20:50,40:70) = NaN;
znan(30:90,5:10) = NaN;
znan(70:75,40:90) = NaN;

tic,z = inpaint_nans(znan,3);toc

tic
k = isnan(znan);
zk = griddata(x(~k),y(~k),z(~k),x(k),y(k));
zg = znan;
zg(k) = zk;
toc

figure
surf(z0)
title 'Original surface'

figure
surf(znan)
title 'Artifacts (large holes) in surface'

figure
surf(zg)
title(['Griddata inpainting (',num2str(sum(isnan(zg(:)))),' NaNs remain)'])

figure
surf(z)
title 'Inpainted surface'

figure
surf(zg-z0)
title 'Griddata error surface'

figure
surf(z-z0)
title 'Inpainting error surface (Note z-axis scale)'

%% Comparison of methods

[x,y] = meshgrid(-1:.02:1);
r = sqrt(x.^2 + y.^2);
z = exp(-(x.^2+ y.^2));

z(r>=0.9) = NaN;

z((r<=.5) & (x<0)) = NaN;

figure
pcolor(z);
title 'Surface provided to inpaint_nans'

% Method 0
tic,z0 = inpaint_nans(z,0);toc

% Method 1
tic,z1 = inpaint_nans(z,1);toc

% Method 2
tic,z2 = inpaint_nans(z,2);toc

% Method 3
tic,z3 = inpaint_nans(z,3);toc

% Method 4
tic,z4 = inpaint_nans(z,4);toc

% Method 5
tic,z5 = inpaint_nans(z,5);toc

figure
surf(z0)
colormap copper
hold on
h = surf(z);
set(h,'facecolor','r')
hold off
title 'Method 0 (Red was provided)'

figure
surf(z1)
hold on
h = surf(z);
set(h,'facecolor','r')
hold off
title 'Method 1 (Red was provided)'

figure
surf(z2)
hold on
h = surf(z);
set(h,'facecolor','r')
hold off
title 'Method 2 (Red was provided) - least accurate, but fastest'

figure
surf(z3)
hold on
h = surf(z);
set(h,'facecolor','r')
hold off
title 'Method 3 (Red was provided) - Slow, but accurate'

figure
surf(z4)
hold on
h = surf(z);
set(h,'facecolor','r')
hold off
title 'Method 4 (Red was provided) - designed for constant extrapolation!'

figure
h = surf(z5);
set(h,'facecolor','y')
hold on
h = surf(z);
set(h,'facecolor','r')
hold off
title 'Method 5 (Red was provided)'


%% 1-d "inpainting" using interp1

x = linspace(0,3*pi,100);
y0 = sin(x);
y = y0;
% Drop out 2/3 of the data
y(1:3:end) = NaN;
y(2:3:end) = NaN;

% inpaint_nans
y_inpaint = inpaint_nans(y,1);

% interpolate using interp1
k = isnan(y);
y_interp1c = y;
y_interp1s = y;
y_interp1c(k) = interp1(x(~k),y(~k),x(k),'cubic');
y_interp1s(k) = interp1(x(~k),y(~k),x(k),'spline');

figure
plot(x,y,'ro',x,y_inpaint,'b+')
legend('sin(x), missing 2/3 points','inpaint-nans','Location','North')

figure
plot(x,y0-y_inpaint,'r-',x,y0-y_interp1c,'b--',x,y0-y_interp1s,'g--')
title 'Inpainting residuals'
legend('Inpaint-nans','Pchip','Spline','Location','North')
