function [proj_lg, BHCalib] = BHCorrection(datafolder, geo, ScanXML, proj_lg,gpuids)
% Entry Function For BH Correction
%   Detailed explanation goes here


disp('BH correction is very memory intensive. The effect is to be justified. ');
typein = input('Continue BH correction? Y or N (recommended): ', 's');
if(contains(typein, 'n', 'IgnoreCase', true))
    BHCalib = NaN;
    disp('BH correction is skipped.');
    return;
end

disp('Beam Hardening Correction is on-going: be patient... ');

% Key calibration information
BHCalib = BHCalibFromXML(datafolder, ScanXML);

% Precompute filter attenuated spectrum: debug pass
BHCalib = BH_SpectrumFilter(BHCalib);

% Precompute bowtie attenuated spectra LUT
BHCalib = BH_SpectrumBowtieLUT(geo, BHCalib);

% Build reference object (water) attanuation LUT
BHCalib = BH_ObjectCalibLUT(BHCalib);

% BH correction via reference object (water)
BHCalib = BH_RemappingFunc(BHCalib);

proj_lg = BH_ObjectRemapping(BHCalib, proj_lg, gpuids);

disp('BH correction is done.')

end

