function [ img ] = Atb( projections,geo,angles,varargin )
%ATB CUDA backprojection operator
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD.
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/blob/master/LICENSE
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
% Lets make 100% sure that data is correct
%% OPtionals

ptype='FDK';
expectedProjectionTypes = {'FDK','matched'};
acceptableOptionName = {'gpuids'};

if nargin > 3
    if any(strcmp(varargin{1}, expectedProjectionTypes))
        ptype = varargin{1};
        [gpuids] = parse_inputs(varargin{2:length(varargin)});
        %[ptype, gpuids] = parse_inputs1(varargin{1}, expectedProjectionTypes, varargin{2:length(varargin)});
    elseif any(strcmp(varargin{1}, acceptableOptionName))
        [gpuids] = parse_inputs(varargin{:});
        %[ptype, gpuids] = parse_inputs1(ptype, expectedProjectionTypes, varargin{:});
    else
        assert(false,'TIGRE:Atb:InvalidInput','Projection type not understood (4th input).');
    end
else
    gpuids = GpuIds();
end


%% image
% gpuArray input: runs on the GPU that holds it (MATLAB's current GPU) and returns a gpuArray.
gpuInput=isa(projections,'gpuArray');
gpuOnHost=false; % run a gpuArray input through the usual host code instead
if gpuInput
    projectionsClass=classUnderlying(projections);
    currentGpu=gpuDevice();
    % As in Python: if 'gpuids' names other GPUs, or several, use them on a host copy.
    gpuOnHost=any(strcmpi(varargin(cellfun(@ischar,varargin)),'gpuids')) && ~isequal(gpuids.devices,int32(currentGpu.Index-1));
else
    projectionsClass=class(projections);
end
assert(strcmp(projectionsClass,'single'),'TIGRE:Atb:InvalidInput','Image should be single type');
assert(isreal(projections),'TIGRE:Atb:InvalidInput','Image should be real (non-complex)');
assert(size(projections,2)>1,'TIGRE:Atb:InvalidInput', 'Projections should be 2D'); %TODO: needed? 
assert(size(projections,3)==size(angles,2),'TIGRE:Atb:InvalidInput', 'Number of projections should match number of angles.'); 
%% Angles
assert(isreal(angles),'TIGRE:Atb:InvalidInput','Angles should be real (non-complex)');
assert(size(angles,1)==1 | size(angles,1)==3 ,'TIGRE:Atb:InvalidInput','Angles should be of size 1xN or 3xN');
angles=double(angles); %in case they were single.
if size(angles,1)==1
   angles=repmat(angles,[3 1]);
   angles(2,:)=0;
   angles(3,:)=0;
end
%% geometry
geo=checkGeo(geo,angles);
assert(isequal([size(projections,2) size(projections,1)],geo.nDetector.'),'TIGRE:checkGeo:BadGeometry','nDetector does not match with provided image size');

%% Thats it, lets call the mex function

if gpuInput && ~gpuOnHost
    try
        img=Atb_mex(projections,geo,angles,ptype, gpuids.devices);
    catch err
        if ~strcmp(err.identifier,'CBCT:MEX:Atb:gpuArray'), rethrow(err); end
        warning('TIGRE:Atb:gpuArray','TIGRE was compiled without gpuArray support (Parallel Computing Toolbox not found when compiling): using a host copy. Recompile to avoid the copies.');
        gpuOnHost=true;
    end
end
if gpuInput && gpuOnHost
    img=gpuArray(Atb_mex(gather(projections),geo,angles,ptype, gpuids.devices));
elseif ~gpuInput
    img=Atb_mex(projections,geo,angles,ptype, gpuids.devices);
end

end

function [gpuids]=parse_inputs(varargin)
    %fprintf('parse_inputs0(varargin (%d))\n', length(varargin));
    if isempty(varargin)
        gpuids = GpuIds();
    else
        % create input parser
        p=inputParser;
        % add optional parameters
        addParameter(p,'gpuids', GpuIds());
        %execute
        parse(p,varargin{:});
        %extract
        gpuids=p.Results.gpuids;
    end
end

