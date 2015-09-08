addpath('bin');

ParamSetting;

%% Recon case 1 - Analytic reconstruction: filtered backprojection
% filter='ram-lak','shepp-logan','cosine', 'hamming', 'hann' : (ramp + additional filter)
param.filter='hann'; 
load proj.mat

proj_filtered = filtering(proj,param);
Reconimg = CTbackprojection(proj_filtered, param);

for i=1:param.nz
    figure(2); imagesc(max(Reconimg(:,:,i),0)); axis off; axis equal; colormap gray; colorbar;
    title(num2str(i));
    pause(0.01);
end