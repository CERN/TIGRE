function Result = ROIMask(img, dia)
%% Cylindrical Mask to Input image matrix
% Date: 2015-10-22
% Author: Yi Du (yi.du@hotmail.com)

[Length, Width, Height] = size(img);

if(Length~=Width)
    error('Length is not equal to Width');
end
Length = Width;

if(Width<dia)
    error('Width is larger than ROI diameter');
end

% diameter -> radius
Rad = dia*0.5;

[XX, YY] = meshgrid((1:Width)-(Width+1)/2, (1:Width)-(Width+1)/2);
XX= XX.^2;
YY = YY.^2;
R2 = Rad^2;

Mask = (XX+YY)-R2;
Mask(Mask>0)=0;
Mask(Mask<0)=1;

MaskMat = repmat(Mask, [1 1 Height]);

Result = img.*MaskMat;

end

