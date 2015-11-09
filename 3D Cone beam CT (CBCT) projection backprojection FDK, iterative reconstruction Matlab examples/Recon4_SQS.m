addpath('bin');

ParamSetting;

%% Recon case 4 -  Iterative reconstruction: SQS
load proj.mat
param.filter = 'hann';

% Initial image is the FDK image
proj_filtered = filtering(proj,param);
init_img = CTbackprojection(proj_filtered, param);

img = init_img;

Norimg = CTbackprojection(proj_exp.*CTprojection(ones(param.nx,param.ny,param.nz,'single'),param), param);

% if you have a scatter sinogram
proj_scatter = 0;

for iter = 1:50
    tic;
    proj_tmp = I0*exp(-CTprojection(img,param));
    proj_diff = (proj_exp./(proj_tmp + proj_scatter) - 1).*proj_tmp;
    
    img_diff = CTbackprojection(proj_diff, param)./Norimg;
    img_diff(isnan(img_diff)) = 0;
    img_diff(isinf(img_diff)) = 0;
    
    img = max(img-img_diff,0);
    toc;
    figure(4); imagesc(max(img(:,:,round(end/2)),0)); axis off; axis equal; colormap gray; colorbar;
    title(['Iteration - ',num2str(iter)]);
    pause(0.1);
end
