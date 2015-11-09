function [ img ] = CTbackprojection( proj, param )
%CTBACKPROJECTION Summary of this function goes here
%   Detailed explanation goes here

img = zeros(param.nx, param.ny, param.nz, 'single');

for i = 1:param.nProj    
    img = img + backprojection(proj(:,:,i),param,i);
end


end

