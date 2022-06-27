function [proj, angles, airnorm] = ProjLoader(datafolder,varargin)
% [proj, angle, blk, projinfo] = BatchReadXim(foldername, varargin)
% Varian .xim files exported from TrueBeam on-board imager/EPID, Halcyon imager
% .xim files are exported under service mode, as a single .zip file
% Once unzipped, all relevant files are organized under the same root folder
% Input:
%       datafolder    : root folder of the upzipped files,
%                       or the subfolder Acquisitions where all xim files are stored
%       thd           : 0(default), no motion lag correcion
%                       thd,  motion correcion threadhold
% Output:
%       proj          : .xim pixel_images frame by frame
%       angles        : KVSourceRtn at each frame
%       projinfo      : proj_*.xim info cell arrays
% Date: 2021-03-23
% Author: Yi Du, yi.du@hotmail.com

thd = 0;
if(nargin == 2)
    thd = varargin{1};
end

fprintf("Dataset in: %s \n", datafolder);

ximfilelist = dir([datafolder filesep '**' filesep 'Proj_*.xim']);

proj = [];
angles = [];
airnorm = [];

%
proj_no = length(ximfilelist);
%
count = 1;
for ii = 1:proj_no
	ximfilename = fullfile(ximfilelist(ii).folder, ximfilelist(ii).name);
    if(~mod(ii,50))
        fprintf("Loading: %d / %d \n", ii, proj_no);
    end
	[page, rtn] = mexReadXim(ximfilename);
	if(~isempty(page))
        % load first page
        if(count==1)
            projinfo{count} = ReadXim(ximfilename, 0);
    		angles(count) = rtn;
        	proj(:,:,count) = rot90(single(page), -1);
            
            airnorm(count) = single(projinfo{count}.properties.KVNormChamber);
            count = count + 1;
        else
            % proj info
            projinfo{count} = ReadXim(ximfilename, 0);
            
            % whether threshold is applied
            if(thd)
                if(abs(rtn-angles(end))>thd)
                    angles(count) = rtn;
                	proj(:,:,count) = rot90(single(page), -1);
                    airnorm(count) = single(projinfo{count}.properties.KVNormChamber);
                    count = count + 1;
                else
                    continue;
                end
            else
                angles(count) = rtn;
            	proj(:,:,count) = rot90(single(page), -1);
                airnorm(count) = single(projinfo{count}.properties.KVNormChamber);
                count = count + 1;
            end
        end
	end
end
fprintf("Loading complete \n");

end
