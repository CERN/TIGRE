function [ imgO ] = imDenoise3D( img,type,varargin )
%IMDENOISE3D removes noise of image with differentmethods
%   tv

if nargin==4
   iter=varargin{2};
   hyper=varargin{1};
else
    iter=50;
    hyper=15;
end

immin=min(img(:));
img=(img-immin);
immax=prctile(img(:),99);
img=img./immax;
if strcmp(type,'TV');  
    imgO=double(permute(tvDenoise(permute(img,[3 2 1]),hyper,iter),[1 2 3]));
end
clear img;
imgO=imgO*immax;
imgO=imgO+immin;
end

