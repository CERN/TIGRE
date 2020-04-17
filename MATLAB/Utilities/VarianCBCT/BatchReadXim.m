function [proj, angles, blk, projinfo] = BatchReadXim(datafolder)
% [proj, angle, blk, projinfo] = BatchReadXim(foldername, varargin)
% Varian .xim files exported from TrueBeam on-board imager/EPID, Halcyon imager
% .xim files are exported under service mode, as a single .zip file
% Once unzipped, all relevant files are organized under the same root folder
% Input:
%       datafolder    : root folder of the upzipped files,
%                       or the subfolder Acquisitions where all xim files are stored
% Output:
%       proj          : .xim pixel_images frame by frame
%       angles         : KVSourceRtn at each frame
%       projinfo      : proj_*.xim info cell arrays
% Date: 2020-04-11
% Author: Yi Du, yi.du@hotmail.com

read_xim_info = 0;
if(nargout==4)
    read_xim_info = 1;
end

fprintf("Dataset in: %s \n", datafolder);

ximfilelist = dir([datafolder filesep '**' filesep 'Proj_*.xim']);

proj = [];
angles = [];

proj_no = length(ximfilelist);

for ii = 1:proj_no
	ximfilename = fullfile(ximfilelist(ii).folder, ximfilelist(ii).name);
    if(~mod(ii,50))
        fprintf("Loading: %d / %d \n", ii, proj_no);
    end
	[page, rtn] = mexReadXim(ximfilename);
	if(~isempty(page))
		proj(:,:,ii) = double(page');
		angles(ii) = rtn;
        if(ii==1)
            projinfo{ii} = ReadXim(ximfilename, 0);
        elseif(read_xim_info)
            projinfo{ii} = ReadXim(ximfilename, 0);
        end
	end
end
fprintf("Loading complete \n");

%% Load Blk
tag_bowtie = strfind(projinfo{1}.properties.KVCollimatorShape, 'None');
% wi Bowtie1
if(isempty(tag_bowtie))
    blkfilestr = dir([datafolder filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'FilterBowtie.xim']);
else
    blkfilestr = dir([datafolder filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'Filter.xim']);
end

blkfilename = fullfile(blkfilestr.folder, blkfilestr.name);
blk = mexReadXim(blkfilename);
blk = double(blk');

%% Don't know why mAs correction doesn't work
%{
% mAs for projection sampling
mAsproj = projinfo{1}.properties.KVMilliAmperes * projinfo{1}.properties.KVMilliSeconds;
% mAs for blk field
blkinfo = ReadXim(blkfilename, 0);
mAsblk = blkinfo.properties.KVMilliAmperes * blkinfo.properties.KVMilliSeconds;
% mAs output correction
gain = mAsproj/mAsblk;
%}
% Emperical value
gain = 2.5; 
blk = blk * gain;
% denoising blk
blk = medfilt2(blk, [5 5]);

end
