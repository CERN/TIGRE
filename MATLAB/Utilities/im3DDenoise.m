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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
if nargin >= 5
    [gpuids] = parse_inputs(varargin{3:length(varargin)});
else
    gpuids = GpuIds();
end
if nargin>=4
    iter=varargin{1};
    hyper=varargin{2};
else
    iter=50;
    hyper=15;
end

if strcmp(type,'TV')
    immin=min(img(:));
    img=(img-immin);
%     immax=prctile(img(:),99);
    immax=max(img(:));

    img=img./(immax+2*eps);
    if ndims(img)==2
        imgO=tvDenoise(cat(3,img,img),hyper,iter,gpuids.devices);
        imgO=imgO(:,:,1);
    else
        imgO=tvDenoise(img,hyper,iter,gpuids.devices);
    end
    clear img;
    
    imgO=imgO*immax;
    imgO=imgO+immin;
else
    error(['Unknown denoising method:  '  type ]);
end
clear img;

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


