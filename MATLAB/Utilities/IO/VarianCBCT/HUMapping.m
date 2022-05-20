function HU3D = HUMapping(img3D, Protocol)
% Pixel value to HU value Mapping
% Method: using CATPhan 604 calibration data
% Input:
%       img3D: reconstructed image matrix in pixel value
%       Protocol: CBCT scan protocol
% Output:
%       HU3D: image matrix in CT HU value
% Date: 2021-09-03
% Author: Yi Du, yi.du@hotmail.com

HU3D = zeros(numel(img3D), 1);
load('D:\MATLABWorkplace\Ongoing_Projects\Fronc_VarianCBCT\demo_Pixel2HU.mat', 'Pixel2HU');

tmp = getfield(Pixel2HU, Protocol);

[xData, yData] = prepareCurveData(tmp(:,1), tmp(:,2));

% 2nd-order Polynominal Fitting :Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft);

HU3D = reshape(fitresult(img3D), size(img3D));

end
