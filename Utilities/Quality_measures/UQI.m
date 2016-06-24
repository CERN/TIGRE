function [uqi]=UQI(real,res)
%Calculate universal quality index (UQI)to evaluate the degree of
% similarity between the reconstructed and phantom images for chosen ROIs.
% Its value ranges from zero to one. A UQI value closer to one suggests
% better similarity to true image.
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

%% Refferences
% UQI in X-rays:
%Ref: Few-view cone-beam CT reconstruction with deformed prior image
%doi: 10.1118/1.4901265

% Original UQI:
%Ref: A universal image quality index
%doi: 10.1109/97.995823

%% Code
%real = exact phantom
%res = reconstructed image
real=real(:);
res=res(:);

N=length(real);

%Mean
meanreal=mean(real);
meanres=mean(res);


%Variance
% varreal=sum((real-meanreal)^2)/(N-1);
varreal=var(real);
varres=var(res);


%Covariance
cova=sum((res-meanres).*(real-meanreal))/(N-1);
% cova=cov(real,res);

front= (2*cova)/(varres+varreal);
back= (2*meanres*meanreal)/(meanres^2+meanreal^2);

% uni= ((2*cova)/((varres^2)+(varreal^2)))*((2*meanres*meanreal)/(((meanres^2)+(meanreal^2))));
uqi=front*back;





end