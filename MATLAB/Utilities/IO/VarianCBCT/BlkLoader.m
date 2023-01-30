function [Blk, Sec, BlkAirNorm] = BlkLoader(datafolder)
% Varian Blk Projection Loader:
% Load all dataset that are needed for reconstruction
% Tested on TrueBeam 2.5 and 2.7
% Date: 2020-04-16
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = '~/your_data_path/varian/2020-01-01_123456/';

%% Scan parameter
[~, ScanXML] = GeometryFromXML(datafolder);

%% True for CC, and False for CW
if(str2double(ScanXML.Acquisitions.StopAngle.Text) - str2double(ScanXML.Acquisitions.StartAngle.Text)>0)
    RtnDirection = 'CC';
else
    RtnDirection = 'CW';
end

%% TrueBeam Version
if( strfind(ScanXML.Attributes.Version, '2.7'))
    Version = 2.7;
elseif( strfind(ScanXML.Attributes.Version, '2.0'))
    Version = 2.0;
else
    error('Error in Version');
end

%% Bowtie
tag_bowtie = strfind(ScanXML.Acquisitions.Bowtie.Text, 'None');

%% Blk File Struct
% wi Bowtie1
if(isempty(tag_bowtie))
    % TB 2.0
    if(Version == 2.0)
        blkfilestr = dir([datafolder filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'FilterBowtie.xim']);
    elseif(Version == 2.7)
        blkfilestr = dir([datafolder filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'FilterBowtie_'  RtnDirection '*.xim']);
    end
else
    blkfilestr = dir([datafolder filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'Filter.xim']);
end

%% Load Blk and AirNorm for TB2.5
if(Version == 2.0)
    filename = fullfile(blkfilestr.folder, blkfilestr.name);
    tmp = mexReadXim(filename);
    Blk = rot90(single(tmp), -1);
    tmp = ReadXim(filename, 0);
    BlkAirNorm = single(tmp.properties.KVNormChamber);
    Sec = [];
    return;
end

%% Cresent Artifacts correction
if(Version == 2.7)
    for ii = 1:length(blkfilestr)
        filename = fullfile(blkfilestr(ii).folder, blkfilestr(ii).name);
        tmp = mexReadXim(filename);
        Blk(:,:,ii) = rot90(single(tmp), -1);
        tmp = ReadXim(filename, 0);
        BlkAirNorm(ii) = single(tmp.properties.KVNormChamber);
        % GantryRtn
        Sec(ii) = single(tmp.properties.GantryRtn);
    end
end

end
