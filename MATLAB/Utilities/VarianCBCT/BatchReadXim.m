function [proj, angle, ximstr] = BatchReadXim(foldername)
% [proj, angle, ximstr] = BatchReadXim(foldername, varargin)
% Varian .xim files exported from TrueBeam on-board imager/EPID, Halcyon imager
% .xim files are exported under service mode, as a single .zip file
% Once unzipped, all relevant files are organized under the same root folder
% Input:
%       foldername    : root folder of the upzipped files,
%                       or the subfolder Acquisitions where all xim files are stored
% Output:
%       proj          : .xim pixel_images frame by frame
%       angles        : GrantryRtn at each frame
%       ximstr        : .xim info cell arrays
% Date: 2020-04-11
% Author: Yi Du, yi.du@hotmail.com

read_xim_info = 0;
if(nargout==3)
    read_xim_info = 1;
end

ximfilelist = dir([foldername filesep '**' filesep 'Proj_*.xim']);

proj = [];
angle = [];
for ii = 1:length(ximfilelist)
	ximfilename = fullfile(ximfilelist(ii).folder, ximfilelist(ii).name);
	[page, rtn] = mexReadXim(ximfilename);
	if(~isempty(page))
		proj(:,:,ii) = page;
		angle(ii) = rtn;
        if(read_xim_info)
            ximstr{ii} = ReadXim(ximfilename, 0);
        end
	end
end

end

