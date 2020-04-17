function [proj, angles, geo] = EasyDataLoader(datafolder)
% Load all dataset that are needed for reconstruction
% Date: 2020-04-16
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = 'E:\BigData\Edge\CBCT_Export\2020-01-09_144244';

%% Load proj and angle
[proj, angles, blk] = BatchReadXim(datafolder);

proj = log(repmat(blk, [1 1 size(proj,3)])./proj);

% Mediate filtering along colume-orth
for ii = 1:size(proj,3)
    proj(:,:,ii) = ordfilt2(proj(:,:,ii), 5, ones(1,9));
end

% in case of abnormlies
proj(isnan(proj)) = 0;
proj(isinf(proj)) = 0;

% all negative to zeros
proj(proj<0) = 0;

% double to single
proj = single(proj);

% degree to rad
angles = angles/180*pi;

% -------------------- to test ------------------
% Gantry Rotation correction: limitation of current FDK
if(angles(end) - angles(1)>0)
    proj = flip(proj, 3);
end

%% Load geometry
geo = GeometryFromXML(datafolder);

end
