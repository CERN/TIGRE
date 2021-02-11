function img=minimizeTV_prior(img,prior,dtvg,ng,ratio,varargin)
%MINIMIZETV_prior MATLAB wrapper for the CUDA stepest descend minimization of TV
% norm with anatomical priors. Note that this does not minimize the TV noise, 
% using the ROF model, this minimizes the TV alone, using gradient descend. 
% Infinite iterations of this code will lead to a flat image or to the prior image.
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
    [gpuids] = parse_inputs(varargin{6:length(varargin)});
else
    gpuids = GpuIds();
end

img=minPICCS(img,prior,dtvg,ng,ratio);

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
