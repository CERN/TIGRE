function [ScanXML,ReconXML] = ReadXimXML(Ximfolder)
% Dependence: xml2struct.m
% Date: 2020-03-28

% Scan.xml
filestr = dir([Ximfolder filesep 'Scan.xml']);

ScanXML = getfield(xml2struct(fullfile(filestr.folder, filestr.name)),...
    'Scan');

% Reconstruction.xml
filestr = dir([Ximfolder filesep '**' filesep 'Reconstruction.xml']);

ReconXML = getfield(xml2struct(fullfile(filestr.folder, filestr.name)),...
    'Reconstruction');

end