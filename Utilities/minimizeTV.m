function img=minimizeTV(img,varargin)
%MINIMIZETV MATLAB wrapper for the CUDA stepest descend minimization of TV
% norm. Note that this does not minimize the TV noise, using the ROF mdoel,
% this minimizes the TV alone. Infinite iterations of this code will lead
% to a flat image.
if nargin==1
    dtdv=1;
    ng=30;
else
    if nargin == 3
        dtvg=varargin{1};
        ng=varargin{2};
    else
        error('Wrogn amount of inputs');
    end
end
img=minTV(img,dtvg,ng);
    img=(permute(img,[3 2 1]));
end