addpath('bin');

ParamSetting;

%% Recon case 3 -  Iterative reconstruction: SART
load proj.mat
param.filter = 'none';

img = zeros(param.nx, param.ny, param.nz,'single');

Norimg = CTbackprojection(CTprojection(ones(param.nx,param.ny,param.nz,'single'),param), param);

for iter = 1:50
    tic;
    proj_diff = CTprojection(img,param)-proj;
    
    img_diff = CTbackprojection(proj_diff, param)./Norimg;
    img_diff(isnan(img_diff)) = 0;
    img_diff(isinf(img_diff)) = 0;
    
    img = max(img-img_diff,0);
    toc;
    figure(4); imagesc(max(img(:,:,round(end/2)),0)); axis off; axis equal; colormap gray; colorbar;
    title(['Iteration - ',num2str(iter)]);
    pause(0.1);
end
