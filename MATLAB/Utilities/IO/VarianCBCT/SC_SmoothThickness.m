function thickness = SC_SmoothThickness(thickness, sccalib, step_du, step_dv)
%% Gaussian Filter to Smooth EstimatedThickness 
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:
%               thickness: Estimated object thickness, i.e., tau(x,y) in Reference
% Output:
%               sccalib: Scatter Calibration Structure
% Date: 2021-05-05
% Author: Yi Du (yi.du@hotmail.com)

% unit: mm
sigma_u = str2double(sccalib.CalibrationResults.Globals.AsymPertSigmaMMu.Text);
sigma_v = str2double(sccalib.CalibrationResults.Globals.AsymPertSigmaMMv.Text);

%% Gaussian filtering
% only one page
if(size(thickness, 3) == 1)
    %% -------------- to check later the coordinate orientation
    thickness = imgaussfilt(thickness, [sigma_v/step_dv sigma_u/step_du]);
% 3D thickness matrix
elseif(size(thickness, 3) >1)
    for ii = 1: size(thickness, 3)
        %% -------------- to check later the coordinate orientation
        thickness(:,:,ii) = imgaussfilt(thickness(:,:,ii), [sigma_v/step_dv sigma_u/step_du]);
    end
end

end
