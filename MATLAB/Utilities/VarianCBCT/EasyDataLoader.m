function [proj, angle, geo] = EasyDataLoader(datafolder)
% EASYDATALOADER Summary of this function goes here
% Detailed explanation goes here

% datafolder = 'E:\BigData\Edge\CBCT_Export\2019-09-03_114227';
%% Load proj and angle
[proj, angle, blk] = BatchReadXim(datafolder);

for ii = 1:length(angle)
    proj(:,:,ii) = log(blk./(proj(:,:,ii) + eps));
end

% remove possible abnormlies
proj(isnan(proj)) = 0;
proj(isinf(proj)) = 0;
proj(isinf(proj)) = 0;

% all negative to zeros
proj(proj<0) = 0;

% double to single
proj = single(proj);

% degree to rad
angle = angle/180*pi;

% -------------------- to test ------------------
% angle starts from zero and increments positively
% limitation of current FDK
if(angle(end) - angle(1)>0)
    angle = angle - angle(1);
else
    angle = flip(angle);
    angle = angle - angle(1);
    proj = flip(proj,3);
end

%% Load geometry
geo = GeometryFromXML(datafolder);

end
