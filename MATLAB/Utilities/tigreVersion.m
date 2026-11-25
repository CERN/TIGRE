function v = tigreVersion()
%TIGREVERSION TIGRE's version, and the versions of what it runs on.
%
%   v = tigreVersion() returns TIGRE's version as text, for example '3.3.1'.
%   tigreVersion() with no output prints the versions of TIGRE, MATLAB and
%   the operating system, and the names of the GPUs. Paste that into bug
%   reports: https://github.com/CERN/TIGRE/issues
%
%   The version is read from pyproject.toml at the root of the TIGRE
%   folder, which MATLAB and Python share.
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
%--------------------------------------------------------------------------
tigre = 'unknown';
here = fileparts(mfilename('fullpath'));            % .../MATLAB/Utilities
toml = fullfile(here, '..', '..', 'pyproject.toml');
if exist(toml, 'file') == 2
    found = regexp(fileread(toml), '^version\s*=\s*"([^"]+)"', 'tokens', 'once', 'lineanchors');
    if ~isempty(found)
        tigre = found{1};
    end
end
if nargout > 0
    v = tigre;
    return
end
fprintf('TIGRE    %s\n', tigre);
fprintf('MATLAB   %s\n', version);
fprintf('OS       %s\n', computer);
try
    names = getGpuNames();
    if isempty(names) || isempty(names{1})
        fprintf('GPUs     none found\n');
    else
        fprintf('GPUs     %s\n', strjoin(names, ', '));
    end
catch err
    fprintf('GPUs     could not query them (%s)\n', err.message);
end
end
