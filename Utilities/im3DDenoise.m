function [ imgO ] = im3DDenoise( img,type,varargin )
%IMDENOISE3D removes noise of image with different methods
%   Currentyl only TV is supported. INput arguments are the iamge, the type
%   of denoising ('TV' only now) and the parameters for the denoising,
%   being number of iterations and hyperparameter currently available. 
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
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
if nargin==4
   iter=varargin{1};
   hyper=varargin{2};
else
    iter=50;
    hyper=15;
end

if strcmp(type,'TV');
    immin=min(img(:));
    img=(img-immin);
    immax=prctile(img(:),99);
    img=img./immax;

    imgO=(permute(tvDenoise(permute(img,[3 2 1]),hyper,iter),[1 2 3]));
    clear img;
    
    imgO=imgO*immax;
    imgO=imgO+immin;
else
    error(['Unknown denoising method:  '  type ]);
end
clear img;

end

