function [rootm]=RMSE(real,res)
%real = exact phantom
%res = reconstructed image
real=real(:);
res=res(:);

N=length(real);
%N = number of voxels
diff=real-res;

rootm=sqrt(sum(diff.^2)/N);



end