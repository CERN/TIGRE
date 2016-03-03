function [cor]=CC(real,res)
%Compute the Pearson correlation coefficient to measure the linear
%dependence between two images


%real = exact phantom
%res = reconstructed image
real=real(:);
res=res(:);

N=length(real);

%compute the mean pixel values of the two images
meanreal=mean(real);
meanres=mean(res);

diffreal=real-meanreal;
diffres=res-meanres;

a=sqrt(sum(diffreal.^2));
b=sqrt(sum(diffres.^2));

cor=sum(diffreal.*diffres)/(a.*b);


end