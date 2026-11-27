function [proj_lg, BHCalib] = BHCorrection(datafolder, geo, ScanXML, proj_lg,gpuids)
% Entry Function For BH Correction
%   Detailed explanation goes here


disp('Beam Hardening Correction is on-going (slow, and only worth it if the images show cupping): be patient... ');

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

