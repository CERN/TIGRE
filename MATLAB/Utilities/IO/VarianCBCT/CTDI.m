function CTDI = CTDI(datafolder)
% Calculate CT Dosimetry Index (cGy) of a CBCT scan 
% CTDI = DoseFactor * (mAs/100) * Num_Projection 
% Reference: TrueBeam Technical Reference Guide II: Imaging
% Date: 2021-11-23
% Author: Yi Du
% Email: yi.du@hotmail.com

% Load ScanXML
if(~exist('ScanXML', 'var'))
    % Read ScanXML
    filestr = dir([datafolder filesep 'Scan.xml']);
    ScanXML = getfield(xml2struct(fullfile(filestr.folder, filestr.name)), 'Scan');
end

% Dose Factor: cGy /100 mAs
DoseFc = str2double(ScanXML.Acquisitions.DoseFactor.Text);
mA = str2double(ScanXML.Acquisitions.Current.Text);
mS = str2double(ScanXML.Acquisitions.PulseLength.Text);
% mAs
mAs = mA*mS/1000;

filedir = dir([datafolder filesep 'Acquisitions' filesep '*' filesep 'Proj_*']);

% FileSize Threshold 10 kB to Exclude Invalid Acquisition 
threshold = 10*1024;
filesize = zeros(length(filedir),1);

for ii=1:length(filedir)
    filesize(ii) = (filedir(ii).bytes >threshold);
end

% Projection Number
nProj = sum(filesize);

% CTDIw, cGy
CTDI = DoseFc * (mAs/100) *nProj;

end

