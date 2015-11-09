addpath('bin');

ParamSetting;

%% Recon case 2 -  Iterative reconstruction: MLEM
load proj.mat
param.filter = 'none';

img = ones(param.nx, param.ny, param.nz,'single');
Norimg = CTbackprojection(ones(param.nu, param.nv, param.nProj, 'single'), param);

for iter = 1:50
    tic;
    proj_ratio = proj./CTprojection(img,param);
    proj_ratio(isnan(proj_ratio)) = 0;
    proj_ratio(isinf(proj_ratio)) = 0;
    
    img_ratio = CTbackprojection(proj_ratio, param)./Norimg;
    img_ratio(isnan(img_ratio)) = 0;
    img_ratio(isinf(img_ratio)) = 0;
    
    img = max(img.*img_ratio,0);
    toc;
    figure(3); imagesc(max(img(:,:,round(end/2)),0)); axis off; axis equal; colormap gray; colorbar;
    title(['Iteration - ',num2str(iter)]);
    pause(0.1);
end
