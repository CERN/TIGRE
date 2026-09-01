function [common, leftover] = parseCommonInputs(proj, geo, angles, argin, algName)
%PARSECOMMONINPUTS Parses the option name/value pairs shared by (almost)
% all TIGRE reconstruction algorithms.
%
% [COMMON, LEFTOVER] = PARSECOMMONINPUTS(PROJ, GEO, ANGLES, ARGIN, ALGNAME)
% parses ARGIN (an algorithm's varargin, as a cell array) for the options
% shared across TIGRE's iterative algorithms -- 'init', 'initimg',
% 'verbose', 'qualmeas', 'gpuids', 'groundtruth' -- and returns them as
% fields of the struct COMMON:
%
%   common.res           - initial image, resolved from 'init'/'initimg'
%   common.verbose        - verbose flag
%   common.QualMeasOpts   - quality measure cell array
%   common.gpuids         - GpuIds object
%   common.gt             - ground truth image (or nan if not given)
%
% Any name/value pairs in ARGIN that are not among the common options
% (e.g. an algorithm-specific option like 'lambda' or 'restart') are
% left untouched, in their original order, in the cell array LEFTOVER,
% so the calling algorithm can parse its own remaining options from
% LEFTOVER exactly as it parsed the full ARGIN before this refactor.
%
% ALGNAME is used only to build 'TIGRE:<ALGNAME>:InvalidInput' error
% identifiers, matching each algorithm's pre-existing error IDs.
%
% NOTE on preserved (not "fixed") behavior: matching every existing
% algorithm's current implementation, if 'init' is set to 'image' but
% 'initimg' is not supplied, COMMON.RES is left as [] without raising an
% error here -- this mirrors the pre-refactor behavior exactly. This is
% arguably a latent bug worth its own separate issue, but changing it is
% out of scope for this structural refactor.
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
% CERN-European Organization for Nuclear Research
% All rights reserved.
%
% License: Open Source under BSD.
% See the full license at
% https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact: tigre.toolbox@gmail.com
% Codes: https://github.com/CERN/TIGRE/
%--------------------------------------------------------------------------

commonOpts = {'init','initimg','verbose','qualmeas','gpuids','groundtruth'};

nVarargs = length(argin);
if mod(nVarargs,2)
    error(['TIGRE:' algName ':InvalidInput'],'Invalid number of inputs')
end

leftover = {};
value = containers.Map('KeyType','char','ValueType','any');
present = containers.Map('KeyType','char','ValueType','logical');

for ii=1:2:nVarargs
    key = lower(argin{ii});
    if ismember(key, commonOpts)
        value(key) = argin{ii+1};
        present(key) = true;
    else
        leftover{end+1} = argin{ii}; %#ok<AGROW>
        leftover{end+1} = argin{ii+1}; %#ok<AGROW>
    end
end

has = @(k) isKey(present, k);

% % % % % % % Verbose
if has('verbose')
    common.verbose = value('verbose');
else
    common.verbose = 1;
end
if ~is2014bOrNewer
    warning(['TIGRE:' algName],'Verbose mode not available for older versions than MATLAB R2014b');
    common.verbose = false;
end

% % % % % % % Init / Initimg
if has('init')
    initVal = value('init');
else
    initVal = 'none';
end
common.res = [];
if (~has('init')) || (ischar(initVal) && strcmp(initVal,'none'))
    common.res = zeros(geo.nVoxel','single');
elseif ischar(initVal) && strcmp(initVal,'FDK')
    common.res = FDK(proj,geo,angles);
elseif ischar(initVal) && strcmp(initVal,'multigrid')
    common.res = init_multigrid(proj,geo,angles);
elseif ischar(initVal) && strcmp(initVal,'image')
    if has('initimg')
        initimgVal = value('initimg');
        if isequal(size(initimgVal),geo.nVoxel')
            common.res = single(initimgVal);
        else
            error(['TIGRE:' algName ':InvalidInput'],'Invalid image for initialization');
        end
    end
    % else: preserved quirk -- res stays [], see NOTE above.
else
    error(['TIGRE:' algName ':InvalidInput'],'Invalid Init option')
end

% % % % % % % Quality measures
if has('qualmeas')
    v = value('qualmeas');
    if iscellstr(v)
        common.QualMeasOpts = v;
    else
        error(['TIGRE:' algName ':InvalidInput'],'Invalid quality measurement parameters');
    end
else
    common.QualMeasOpts = {};
end

% % % % % % % GPU ids
if has('gpuids')
    common.gpuids = value('gpuids');
else
    common.gpuids = GpuIds();
end

% % % % % % % Ground truth
if has('groundtruth')
    common.gt = value('groundtruth');
else
    common.gt = nan;
end

end
