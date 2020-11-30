function [proj,geo, angles ] = VarianDataLoader(datafolder, varargin)
% VarianDataLoader   Loads Varian CBCT projection, geomtry and angles data
%
%   Optional parameter: Motion lag correction. Default True. 
%
% Load all dataset that are needed for reconstruction
% Date: 2020-04-16
% Author: Yi Du (yi.du@hotmail.com)
% datafolder = 'E:\BigData\Edge\CBCT_Export\2020-01-09_144244';

%% Load geometry
[geo, ScanXML] = GeometryFromXML(datafolder);

%% Motion lag correcion
thd = 0;
if(~isempty(varargin)&&(varargin{1}))
    thd = str2double(ScanXML.Acquisitions.Velocity.Text)...
        ./str2double(ScanXML.Acquisitions.FrameRate.Text);
    thd = thd *0.95;
end

%% Load proj and angle
[proj, angles, blk] = BatchReadXim(datafolder, thd);

%% Logarithmic calculation
proj = log(repmat(blk, [1 1 size(proj,3)])./proj);

%% Mediate filtering along colume-orth
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

%% Gantry Rotation correction
% Clockwise
if(angles(end) - angles(1)>0)
    proj = flip(proj, 3);
% Counter-clockwise -> Clockwise
else
    angles = flip(angles);
end

end
