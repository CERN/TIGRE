function gform = SC_FormFunc(sccalib, dugd, dvgd)
%% Thickness-based Multiple Group Form Function kernels
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Date: 2021-05-18
% Author: Yi Du

%% group number
ngroup = length(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel);

mm2cm = 0.1;
dugd = dugd * mm2cm;
dvgd = dvgd * mm2cm;
grid2 = dugd.^2 +dvgd.^2;

%% Kernel form function groups
% Form function groups
gform = [];

for ii=1:ngroup
    tmp = sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{ii}.ObjectScatterFit;
    % unit: cm
    sigma1 = str2double(tmp.sigma1.Text);
    sigma2 = str2double(tmp.sigma2.Text);
    % unitless
    B = str2double(tmp.B.Text);

    % Form Function
    gform(:,:,ii) = exp( -0.5 * grid2 /(sigma1^2)  ) + B * exp( -0.5 * grid2 /(sigma2^2) );
end


end
