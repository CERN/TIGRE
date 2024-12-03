function thickness = SC_ThicknessEstimator(blk, page, sccalib, step_du, step_dv)
%% Estimate Water-Equivalent Thickness (2D)
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:
%               Blk(u,v): Total intensity, i.e., I_0
%               BlkAirNorm: 
%               Prm(u,v): Primary intensity, i.e., I_p
%               AirNorm: Primary intensity airnorm chamber value
%               sccalib: Scatter Calibration Structure
% Output:
%               thickness(u,v): Estimated object thickness, i.e., tau(x,y) in Reference
% Date: 2021-05-05
% Author: Yi Du (yi.du@hotmail.com)

% mu H2O = 0.02 /mm
muH2O = str2double(sccalib.CalibrationResults.Globals.muH2O.Text);

% unit mm
tmp = blk./page;
tmp(tmp<0)=0.0001;
thickness = log(tmp) /muH2O;

% fill holes by interpolation
thickness(thickness<0) = NaN;

thickness = single(inpaint_nans(double(thickness), 2));

%% Smooth the estimated thickness
% thickness(vv, uu, ntheta)
thickness = SC_SmoothThickness(thickness, sccalib, step_du, step_dv);

end
