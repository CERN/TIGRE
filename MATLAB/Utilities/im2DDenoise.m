function img = im2DDenoise(img, type, kernelsize)
%img2DDenoise removes noise of image with 2D kernel-based methods
%   Currentyl only median and average filter is supported.
%       type: 'median' or 'average'
%       kernelsize: filter kernel size
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Date: 2021-09-21
% Author: Yi Du, 
% Contact: yi.du@hotmail.com


if strcmp(type,'median')
    SliceNum = size(img,3);
    for kk =1:SliceNum
        img(:,:,kk) = medfilt2(img(:,:,kk), [kernelsize, kernelsize]);
    end
elseif strcmp(type,'average')
    h = fspecial('average', kernelsize);
    img = imfilter(img, h, 'same');
end    

end

