classdef defaultGeometry < Geometry
%%geo = defaultGeometry() generates a default geometry for tests.
% Optional parameters
%
% 'mode'      : 'parallel' or 'cone'
% 'nVoxel'    : 3x1 matrix of size of the image
% 'nDetector' : 2x1 matrix of size of the detector
%
    methods
        function obj = defaultGeometry(varargin)
            [mode, nVoxel, nDetector]=parse_inputs(varargin{:});
            if strcmp(mode,'cone')
                obj = construct_conebeam_geometry(obj, nVoxel, nDetector);
            else
                obj = construct_parallel_geometry(obj, nDetector);
            end
        end
        
        function obj = construct_conebeam_geometry(obj, nVoxel, nDetector)
            %% Example
            % VARIABLE                                   DESCRIPTION                    UNITS
            %-------------------------------------------------------------------------------------
            % Distances
            obj.DSD = 1536;                             % Distance Source Detector      (mm)
            obj.DSO = 1000;                             % Distance Source Origin        (mm)
            % Detector parameters
            obj.nDetector=nDetector;					% number of pixels              (px)
            obj.sDetector=[512*0.8;512*0.8];            % total size of the detector    (mm)
            obj.dDetector=obj.sDetector./obj.nDetector;					% size of each pixel            (mm)

            % Image parameters
            obj.nVoxel=nVoxel;                          % number of voxels              (vx)
            obj.sVoxel=[256;256;256];                   % total size of the image       (mm)
            obj.dVoxel=obj.sVoxel./obj.nVoxel;          % size of each voxel            (mm)
            % Offsets
            obj.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
            obj.offDetector=[0; 0];                     % Offset of Detector            (mm)

            % Auxiliary
            obj.accuracy=1;
            obj.mode='cone';
        end
        function obj = construct_parallel_geometry(obj, nDetector)
                   %% Example
            % VARIABLE                                   DESCRIPTION                    UNITS
            %-------------------------------------------------------------------------------------
            % Distances
            obj.DSD = nDetector(1)*4;                             % Distance Source Detector      (mm)
            obj.DSO = nDetector(1)*2;                             % Distance Source Origin        (mm)
            % Detector parameters
            obj.nDetector=nDetector;					% number of pixels              (px)
            obj.dDetector=[1; 1]; 					% size of each pixel            (mm)
            obj.sDetector=obj.nDetector; % total size of the detector    (mm)
            % Image parameters
            obj.nVoxel=[nDetector(1);nDetector(1);nDetector(2)];                   % number of voxels              (vx)
            obj.sVoxel=obj.nVoxel;                   % total size of the image       (mm)
            obj.dVoxel=[1;1;1];          % size of each voxel            (mm)
            % Offsets
            obj.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
            obj.offDetector=[0; 0];                     % Offset of Detector            (mm)

            % Auxiliary
            obj.accuracy=1;
            obj.mode='parallel'; 
        end
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
    warning('In Parallel mode nVoxel(2:3) is generally equal to nDetector. Consider setting them equal for better recosntruction quality.');
end

end