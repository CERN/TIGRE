function BHCalib = BH_SpectrumBowtieLUT_deprecated(geo, BHCalib)
%SPECTRAPROCESS Summary of this function goes here
%   Detailed explanation goes here

%% Thin window filter has been applied to the spectrum already

% Source-to-detector distance: mm
SDD = geo.DSD;
% Source-to-bowtie distance: mm
SBD = BHCalib.bowtie.distance;

%% [u,v] vector: mm
offset=geo.offDetector;
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2);

%% to bowtier surface
ub = us * SBD/SDD;
% tranverse length vector: mm
ul = interp1(BHCalib.bowtie.uu, BHCalib.bowtie.thickness, ub);
% tranverse length grid: mm
ulgd = repmat(ul, [length(vs), 1]);

%% bowtie material attenuation look-up table
[Min, Max] = bounds(ulgd(:));
% sampling length
n_samples=100;
sl = linspace(Min, Max, n_samples);
% LUT: [sl, miu_bowtie]
atten_tab = zeros(n_samples, length(BHCalib.bowtie.ac));

% I/I_0 = exp(- length * miu): dependent on incident energy
for ii = 1:n_samples
    atten_tab(ii,:) = exp(-sl(ii).* BHCalib.bowtie.ac');
end

%% bowtie-attenuated spectra look-up table
kVp = length(BHCalib.source.spec);
% specLUT(sampling length, energy bin)
BHCalib.bowtie.specLUT = atten_tab(:,1:kVp).*BHCalib.source.spec;
% miu
BHCalib.object.ac =BHCalib.object.ac(1:kVp);

%% sampling length
BHCalib.bowtie.sl = sl;

%% tranvers length of each detector unit
BHCalib.bowtie.ulgd = ulgd;

end
