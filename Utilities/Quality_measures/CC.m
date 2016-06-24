function [cor]=CC(real,res)
%Compute the Pearson correlation coefficient to measure the linear
% dependence between two images
%
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
% Coded by:           Manasavee Lohvithee
%--------------------------------------------------------------------------

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