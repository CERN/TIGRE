function geo=defaultGeometry(varargin)
%%geo = defaultGeometry() generates a default geomtry for tests.
% Optional parameters
%
% 'mode'      : 'parallel' or 'cone'
% 'nVoxel'    : 3x1 matrix of size of the image
% 'nDetector' : 2x1 matrix of size of the detector
%

[mode,nVoxel,nDetector]=parse_inputs(varargin{:});


if strcmp(mode,'cone')
    %% Example
    % VARIABLE                                   DESCRIPTION                    UNITS
    %-------------------------------------------------------------------------------------
    % Distances
    geo.DSD = 1536;                             % Distance Source Detector      (mm)
    geo.DSO = 1000;                             % Distance Source Origin        (mm)
    % Detector parameters
    geo.nDetector=nDetector;					% number of pixels              (px)
    geo.sDetector=[512*0.8;512*0.8];            % total size of the detector    (mm)
    geo.dDetector=geo.sDetector./geo.nDetector;					% size of each pixel            (mm)

    % Image parameters
    geo.nVoxel=nVoxel;                          % number of voxels              (vx)
    geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
    geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
    % Offsets
    geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
    geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
    
    % Auxiliary
    geo.accuracy=1;
    geo.mode=mode;
else 
        %% Example
    % VARIABLE                                   DESCRIPTION                    UNITS
    %-------------------------------------------------------------------------------------
    % Distances
    geo.DSD = nDetector(1)*4;                             % Distance Source Detector      (mm)
    geo.DSO = nDetector(1)*2;                             % Distance Source Origin        (mm)
    % Detector parameters
    geo.nDetector=nDetector;					% number of pixels              (px)
    geo.dDetector=[1; 1]; 					% size of each pixel            (mm)
    geo.sDetector=geo.nDetector; % total size of the detector    (mm)
    % Image parameters
    geo.nVoxel=[nDetector(1);nDetector(1);nDetector(2)];                   % number of voxels              (vx)
    geo.sVoxel=geo.nVoxel;                   % total size of the image       (mm)
    geo.dVoxel=[1;1;1];          % size of each voxel            (mm)
    % Offsets
    geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
    geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
    
    % Auxiliary
    geo.accuracy=1;
    geo.mode=mode;
end
end
function [mode,nVoxel,nDetector]=parse_inputs(varargin)
% create input parser
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

if size(nDetector,1)==1
    nDetector=nDetector.';
end
if size(nVoxel,1)==1
    nVoxel=nVoxel.';
end
if strcmp(mode,'parallel')&&(all(nVoxel(2:3)~= nDetector))
    warning('In Parallel mode nVoxel(2:3) is generally equal to nDetector. Consider setting them equal for better reconstruction quality.');
end

end