function BHCalib = BH_RemappingFunc(BHCalib)
% remaps the non-linear projection signals via linear fitted models
%
% SYNOPSIS: BHCalib = BH_RemappingFunc(BHCalib)
%
% INPUT BH Calibration Structure
%
% OUTPUT BHCalib linear mapping functions
%
% REMARKS
%
% created with MATLAB ver.: 9.8.0.1873465 (R2020a) Update 8 on Microsoft Windows 10 Pro Version 10.0 (Build 18363)
%
% created by: Yi Du
% DATE: 16-May-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bowtie_sl = length(BHCalib.bowtie.sl);
object_sl = length(BHCalib.object.sl);

% control percentage for linear fitting: 5% (default)
cp = 0.2;

fitted_sl = round(object_sl *cp);
object_sl = BHCalib.object.sl(1:fitted_sl);

% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = 0.01;

cof = zeros(bowtie_sl,1);
rsquare = cof;

%% Calculate Fitting Parameters
disp('Calculating Linear Fitting Model: ')
for ii = 1:bowtie_sl    
    [xData, yData] = prepareCurveData(object_sl, BHCalib.object.calibLUT(ii, 1:fitted_sl));
    % Fit model to data.
    [func_ft, gof] = fit( xData, yData, ft, opts );
    cof(ii) = func_ft.a;
    rsquare(ii) = gof.rsquare;
end

%% Fitting Model to Return 
BHCalib.object.linear_ft_cof = cof;

%{
plot(BHCalib.object.sl, BHCalib.object.calibLUT(ii,:), 'r'), grid on; hold on;
plot(BHCalib.object.sl,cof(ii).*BHCalib.object.sl, 'b'), grid on;
plot(BHCalib.object.sl,cof(ii).*BHCalib.object.sl - BHCalib.object.calibLUT(ii,:), 'k'), grid on;
%}

end
