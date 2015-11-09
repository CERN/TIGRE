addpath('bin');

ParamSetting;
param.gpu=0;
%% Make measurement - projection
load img128.mat % Ground-truth image
img=ones(512,512,512);
NoiseOn = 0;    % 0- without noise, 1- with noise

proj = CTprojection(img,param);

if NoiseOn == 1
    I0 = 10000; % higher I0 - small noise, smaller I0 - high noise (You can change this)
    proj_exp = poissrnd( max(I0*exp(-proj),1) );
    proj = -log(min(proj_exp/I0, 1) );
else
    I0 = 1;
    proj_exp = I0*exp(-proj);
end

save proj.mat proj proj_exp I0

% plot
for i=1:param.nProj
    figure(1); imagesc(max(proj(:,:,i)',0)); axis image; axis equal; colormap gray; colorbar;
    title(['Degree : ',num2str(param.deg(i))]);
    pause(0.01);
end
    
