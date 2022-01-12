function img=minimizeAwTV(img,varargin)
%MINIMIZEAWTV MATLAB wrapper for the CUDA steepest descend minimization of 
% adaptive-weighted TV norm. Note that this does not minimize the TV noise, using the ROF mdoel,
% this minimizes the TV alone. Infinite iterations of this code will lead
% to a flat image.
%--------------------------------------------------------------------------
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
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
if nargin >= 6
    [gpuids] = parse_inputs(varargin{4:length(varargin)});
else
    gpuids = GpuIds();
end
if nargin==1
    dtvg=1;
    ng=30;
    delta=-0.005;
else
    if nargin == 6
        dtvg=varargin{1};
        ng=varargin{2};
        delta=varargin{3};
    else
        error('Wrong amount of inputs, 1 or 6 expected');
    end
end

img=AwminTV(img,dtvg,ng,delta,gpuids.devices);

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

