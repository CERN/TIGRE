function img=minimizeTV(img,varargin)
%MINIMIZETV MATLAB wrapper for the CUDA stepest descend minimization of TV
% norm. Note that this does not minimize the TV noise, using the ROF mdoel,
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
if nargin==1
    dtvg=1;
    ng=30;
else
    if nargin == 3
        dtvg=varargin{1};
        ng=varargin{2};
    else
        error('Wrogn amount of inputs');
    end
end
if ndims(img)==3
    img=minTV(img,dtvg,ng);
else
    for ii=1:ng
       gd=gradientTVnorm(img,'backward');
       img=img-dtvg*gd/sqrt(sum(gd(:).^2));
    end
end
end