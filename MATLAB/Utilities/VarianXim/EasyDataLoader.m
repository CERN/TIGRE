function [proj, angle, geo, ximinfo, ScanXML, ReconXML] = EasyDataLoader(datafolder, varargin)
% EASYDATALOADER Summary of this function goes here
% Detailed explanation goes here

% datafolder = 'E:\BigData\Edge\CBCT_Export\2019-09-03_121026';
%% Learn to use inputParser
p=inputParser;
% add optional parameters
addParameter(p,'mode','cone',@(x)(ischar(x)&&(strcmp(x,'parallel')||strcmp(x,'cone'))));
addParameter(p,'nVoxel',   [256;256;256],@(x)(isnumeric(x)&&all(floor(x)==x)&&max(size(x))==3));
addParameter(p,'nDetector',[256;256]    ,@(x)(isnumeric(x)&&all(floor(x)==x)&&max(size(x))==2));

%execute
parse(p,varargin{:});
%extract
mode=p.Results.mode;
nVoxel=p.Results.nVoxel;
nDetector=p.Results.nDetector;

%% Load proj and angle
filestr = dir([datafolder filesep '**' filesep 'Proj_*.xim']);

proj = [];
angle = [];
for ii = 1:length(filestr)
	filename = fullfile(filestr(ii).folder, filestr(ii).name);
	[page, rtn] = mexReadXim(filename);
	if(~isempty(page))
		proj(:,:,ii) = page;
		angle(ii) = rtn;
        ximinfo{ii} = ReadXim(filename, 0);
	end
end

%% Generate geo



end
