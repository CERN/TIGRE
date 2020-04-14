function [proj, angle, blk, projinfo] = BatchReadXim(foldername)
% [proj, angle, blk, projinfo] = BatchReadXim(foldername, varargin)
% Varian .xim files exported from TrueBeam on-board imager/EPID, Halcyon imager
% .xim files are exported under service mode, as a single .zip file
% Once unzipped, all relevant files are organized under the same root folder
% Input:
%       foldername    : root folder of the upzipped files,
%                       or the subfolder Acquisitions where all xim files are stored
% Output:
%       proj          : .xim pixel_images frame by frame
%       angles        : GrantryRtn at each frame
%       projinfo        : proj_*.xim info cell arrays
% Date: 2020-04-11
% Author: Yi Du, yi.du@hotmail.com

read_xim_info = 0;
if(nargout==4)
    read_xim_info = 1;
end

ximfilelist = dir([foldername filesep '**' filesep 'Proj_*.xim']);

proj = [];
angle = [];
for ii = 1:10% length(ximfilelist)
	ximfilename = fullfile(ximfilelist(ii).folder, ximfilelist(ii).name);
	[page, rtn] = mexReadXim(ximfilename);
	if(~isempty(page))
		proj(:,:,ii) = page;
		angle(ii) = rtn;
        if(ii==1)
            projinfo{ii} = ReadXim(ximfilename, 0);
        elseif(read_xim_info)
            projinfo{ii} = ReadXim(ximfilename, 0);
        end
	end
end

%% Load Blk
tag_bowtie = strfind(projinfo{1}.properties.KVCollimatorShape, 'None');
% wi Bowtie1
if(isempty(tag_bowtie))
    blkfilestr = dir([foldername filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'FilterBowtie.xim']);
else
    blkfilestr = dir([foldername filesep 'Calibrations' filesep 'AIR-*' filesep '**' filesep 'Filter.xim']);
end

blkfilename = fullfile(blkfilestr.folder, blkfilestr.name);
blk = mexReadXim(blkfilename);

end
