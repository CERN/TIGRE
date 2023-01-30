function kernel = SC_ASGkernel(sccalib, geo, dus, dvs)
%% Anti-sctter grid response function
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:    
%           geo: geometry structure
%           dus: downsampled u vector
%           dvs: downsampled v vector
% Output:
%           kernel: anti-scatter grid
%
% Date: 2021-05-04
% Author: Yi Du (yi.du@hotmail.com)

%% Geometry
% unit mm
DSD = geo.DSD;

%% Anti-Scatter Grid along X direction
% vs dimension
gamma = abs(rad2deg(atan(dvs'./DSD)));

%% Transmission modelling
k = -0.15;
b = 1;

t_ratio = k.*abs(dvs'/10) + b;

%% Kernel: [nv, nu]
kernel = repmat(t_ratio, [1,length(dus)]);
efficiency = str2double(sccalib.CalibrationResults.ObjectScatterModels.ObjectScatterModel{1}.GridEfficiency.LamellaTransmission.Text);
kernel(kernel < efficiency) = efficiency;

end
