function BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib)
%SPECTRAPROCESS Summary of this function goes here
%   Detailed explanation goes here

%% Thin window filter has been applied to the spectrum already
% Source-to-detector distance: mm
SDD = mean(geo.DSD);
% Source-to-bowtie distance: mm
SBD = BHCalib.bowtie.distance;
% bowtie height = 55 mm, ref to TrueBeam Technical Reference: rays above or below miss it
bowtie_height = 55;

%% [u,v] vector: mm
offset = geo.offDetector(:,1);
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);

%% Fitting to Model Bowtie Profile: 
[xData, yData] = prepareCurveData(BHCalib.bowtie.uu, BHCalib.bowtie.thickness);

%% Fitting Model: Smoothing Spline 
% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft);
disp(['Bowtie Fitting Goodness (R^2): ' num2str(gof.rsquare)])

%% Length of each detector ray inside the bowtie: mm
% The bowtie fills the region between SBD and SBD+thickness(y) from the source, where y is the
% lateral (u) position at that distance, and |z| <= bowtie_height/2. March along each ray
% through it (this replaces voxelising the bowtie and forward projecting it with Ax).
t_max = ceil(max(BHCalib.bowtie.thickness));
n_steps = 4000;
step = t_max/n_steps;
depth = ((1:n_steps)' - 0.5)*step;   % mm past the bowtie front face, along the central ray
dist = SBD + depth;                  % distance from the source
y = dist*us/SDD;                     % lateral position at that distance: [n_steps, nu]
thick = reshape(fitresult(y(:)), size(y));
thick = min(max(thick, 0), t_max);
% length inside the bowtie up to each depth, along the beam axis
cum = [zeros(1, numel(us)); cumsum(double(depth < thick), 1)]*step;
% rays leave the bowtie above/below it: only depths with |z| <= bowtie_height/2 count
dist_max = (bowtie_height/2)*SDD./abs(vs);
dist_max(vs == 0) = Inf;
k = sum(dist(:) <= dist_max(:)', 1);  % number of depth samples inside, per detector row
ulgd = cum(k+1, :);                   % [nv, nu], as Ax returns it
% along the ray itself
ulgd = ulgd .* sqrt(SDD^2 + us.^2 + vs(:).^2)/SDD;


%% bowtie material attenuation look-up table
[Min, Max] = bounds(ulgd, 'all');
% sampling length
sl = linspace(Min, Max, 2000);
% LUT: [sl, miu_bowtie]
atten_tab = zeros(length(sl), length(BHCalib.bowtie.ac));

% I/I_0 = exp(- length * miu): dependent on incident energy
for ii = 1:length(sl)
    atten_tab(ii,:) = exp(-sl(ii).* BHCalib.bowtie.ac');
end

% kVp
kVp = length(BHCalib.filter.spec);
% miu
BHCalib.object.ac =BHCalib.object.ac(1:kVp);
%% bowtie-attenuated spectra look-up table
% specLUT(sampling length, energy bin)
BHCalib.bowtie.specLUT = atten_tab(:,1:kVp).*BHCalib.filter.spec;

%% sampling length vector
BHCalib.bowtie.sl = sl;

%% tranvers length map of each detector unit
BHCalib.bowtie.ulgd = ulgd;

end
