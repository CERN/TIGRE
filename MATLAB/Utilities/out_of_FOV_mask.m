function result = out_of_FOV_mask(img, geo, dia)
% ROI_mask=Cylindrical Mask to Input image matrix
% Date: 2015-10-22
% Author: Yi Du (yi.du@hotmail.com)
% Modified: Ander Biguri

[Length, Width, ~] = size(img);

% diameter -> radius
Rad = dia*0.5;

[XX, YY] = meshgrid((1:Width)-(Width+1)/2-geo.offOrigin(2)./geo.dVoxel(2), (1:Length)-(Length+1)/2-geo.offOrigin(1)./geo.dVoxel(1));
XX= XX.^2;
YY = YY.^2;
R2 = Rad^2;

mask = (XX+YY)-R2;
mask(mask>0)=0;
mask(mask<0)=1;

result = bsxfun(@times,img,mask);

end

