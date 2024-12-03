function [BHCalib, BHCalibXML] = BHCalibFromXML(datafolder, ScanXML)
%% Load Calibration.xml for Beam Hardening Correction
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Input:
%               datafolder: Varian CBCT scan data folder
%               ScanXML: scan geometry structure
% Output:
%               BHCalib: key information for Beam Hardening Calibration
%               BHCalibXML: Beam Hardening Calibration Structure
% Date: 2021-06-01
% Author: Yi Du (yi.du@hotmail.com)

srcfilename = [datafolder filesep 'Calibrations' filesep 'ASC' filesep 'Factory' filesep 'Calibration.xml'];

%% Export as struct
tmp = xml2struct(srcfilename);
BHCalibXML = tmp.Calibration;

%% Restructure BHCalibXML: spectrum data
% Energy Bin Wid
BHCalibXML.BHCalib.energybin = str2double(BHCalibXML.CalibrationResults.EnergyBinWidth.Text);

for ii =1:length(BHCalibXML.CalibrationResults.Spectra.SpectrumProperties)
    % scan voltage
    tf = strcmpi(ScanXML.Acquisitions.Voltage.Text, BHCalibXML.CalibrationResults.Spectra.SpectrumProperties{ii}.Voltage.Text);
    if(tf)
        voltage = str2double(ScanXML.Acquisitions.Voltage.Text);
        tmp = cell2mat(BHCalibXML.CalibrationResults.Spectra.SpectrumProperties{ii}.Flux.float);
        for jj = 1:length(tmp)
            flux(jj) = str2double(tmp(jj).Text);
        end
    % voltage
    BHCalibXML.BHCalib.voltage = voltage;
    % flux from source
    BHCalibXML.BHCalib.source.flux = flux;
    % normalized spectrum from source
    BHCalibXML.BHCalib.source.spec = flux./max(flux(:));
    break;
    end
end

BHCalibXML.BHCalib.source.kV = 1:BHCalibXML.BHCalib.energybin:BHCalibXML.BHCalib.voltage;

%% Restructure BHCalibXML: Filter type and thickness
for ii =1:length(BHCalibXML.CalibrationResults.Filters.FilterProperties)
    % Filter type
    tf = strcmpi(ScanXML.Acquisitions.KVFilter.Text, BHCalibXML.CalibrationResults.Filters.FilterProperties{ii}.Id.Text);
    if(tf)
        % filter name        
        BHCalibXML.BHCalib.filter.name = BHCalibXML.CalibrationResults.Filters.FilterProperties{ii}.Id.Text;
        % filer material
        BHCalibXML.BHCalib.filter.material = BHCalibXML.CalibrationResults.Filters.FilterProperties{ii}.Material.Text;
        % filer thickness (unit: mm)
        BHCalibXML.BHCalib.filter.thickness = str2double(BHCalibXML.CalibrationResults.Filters.FilterProperties{ii}.Thickness.Text);
        break;
    end        
end

%% Restructure BHCalibXML: bowtie coordinates and profiles
% Bowtie type in specific scan is defined in ScanXML
for ii =1:length(BHCalibXML.CalibrationResults.Bowties.BowtieProperties)
    % may be some items are comments rather than structured bowtie info
    if(~isfield(BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}, 'Id'))
        continue;
    end
    tf = strcmpi(ScanXML.Acquisitions.Bowtie.Text, BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}.Id.Text);
    if(tf)
        % Bowtie Id
        BHCalibXML.BHCalib.bowtie.name = BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}.Id.Text;
        % Bowtie distance in mm
        BHCalibXML.BHCalib.bowtie.distance = str2double(BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}.SourceDist.Text);
        if(isfield(BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}, 'Coordinate'))
            % material
            BHCalibXML.BHCalib.bowtie.material = BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}.Material.Text;
            % uu, profile
            tmpuu = cell2mat(BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}.Coordinate.float);
            tmpprofile = cell2mat(BHCalibXML.CalibrationResults.Bowties.BowtieProperties{ii}.Profile.float);
            for jj = 1:length(tmpuu)
                coordinates(jj) = str2double(tmpuu(jj).Text);
                profiles(jj) = str2double(tmpprofile(jj).Text);
            end
            % Coordinates: u-direction, mm
            BHCalibXML.BHCalib.bowtie.uu = coordinates;
            % thickness: mm
            BHCalibXML.BHCalib.bowtie.thickness = profiles;
        end
    end
end

%% Restructure BHCalibXML: object material
% BH calibration material: water in default
BHCalibXML.BHCalib.object.material = BHCalibXML.CalibrationResults.ObjectMaterialId.Text;
% maximal thickness for object calibration: mm
BHCalibXML.BHCalib.object.thicknessmax = str2double(BHCalibXML.CalibrationResults.ObjectThicknessMax.Text);
% thickness sampling delta: mm
BHCalibXML.BHCalib.object.delta = str2double(BHCalibXML.CalibrationResults.ObjectThicknessDelta.Text);

%% Restructure BHCalibXML:  material coefficient database
% attenuation coeeficient: /mm
for ii =1:length(BHCalibXML.CalibrationResults.Materials.MaterialProperties)
    % Filter material
    filter_tf = strcmpi(BHCalibXML.BHCalib.filter.material, BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text);
    if(filter_tf)
        ac = [];
        % attenuation coefficients: /mm
        tmp = cell2mat(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.AttenuationCoefficient.float);
        for jj = 1:length(tmp)
            ac(jj) = str2double(tmp(jj).Text);
        end
        BHCalibXML.BHCalib.filter.ac = ac;
    end
    
    % Object material
    obj_tf = strcmpi(BHCalibXML.BHCalib.object.material, BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text);
    if(obj_tf)
        ac = [];
        % attenuation coefficients: /mm
        tmp = cell2mat(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.AttenuationCoefficient.float);
        for jj = 1:length(tmp)
            ac(jj) = str2double(tmp(jj).Text);
        end
        BHCalibXML.BHCalib.object.ac = ac;
    end
        
    % Bowtie material
    bowtie_tf = strcmpi(BHCalibXML.BHCalib.bowtie.material, BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text);
    if(bowtie_tf)
        ac = [];
        % attenuation coefficients: /mm
        tmp = cell2mat(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.AttenuationCoefficient.float);
        for jj = 1:length(tmp)
            ac(jj) = str2double(tmp(jj).Text);
        end
        BHCalibXML.BHCalib.bowtie.ac = ac;
    end
    
    % Bone material for further correction
    bone_tf = contains(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text, 'bone', 'IgnoreCase',true);
    if(bone_tf)
        ac = [];
        % bone material name
        BHCalibXML.BHCalib.bone.material = BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text;   
        % attenuation coefficients: /mm
        tmp = cell2mat(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.AttenuationCoefficient.float);
        for jj = 1:length(tmp)
            ac(jj) = str2double(tmp(jj).Text);
        end
        BHCalibXML.BHCalib.bone.ac = ac;
    end
    
    % Scitillator material for detection response (PaxScan 4030CB: CsI: Ti)
    scintillator_tf = contains(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text, 'cesium', 'IgnoreCase', true);
    if(scintillator_tf)
        ac = [];
        % detector scintillator material
        BHCalibXML.BHCalib.scintillator.material = BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.Id.Text;   

        % CsI thickness: this parameter is not important at all.
        BHCalibXML.BHCalib.scintillator.thickness = 0.6;   
        
        % attenuation coefficients: /mm
        tmp = cell2mat(BHCalibXML.CalibrationResults.Materials.MaterialProperties{ii}.EnergyAbsorptionCoefficient.float);
        for jj = 1:length(tmp)
            ac(jj) = str2double(tmp(jj).Text);
        end
        % energy absorption coefficient
        BHCalibXML.BHCalib.scintillator.ac = ac;
    end
    
end

%% Key Information for further beam hardening calibration
BHCalib = BHCalibXML.BHCalib;
end
