function [res,errorL2]=FDK(proj,geo,alpha,filter)
%TODO docs FDK
% 
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
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Kyungsang Kim, modified by Ander Biguri 
%--------------------------------------------------------------------------


%Input is data,geosize,angles

if nargin<4
    geo.filter='ram-lak'; 
else
    geo.filter=filter;
end

if size(geo.offDetector,2)==1
    offset=repmat(geo.offDetector,[1 length(alpha)]);
else
    offset=geo.offDetector;
end


%% Weight
%proj=data
proj=permute(proj,[2 1 3]);
for ii=1:length(alpha)
    
    us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1,ii);
    vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2,ii);
    [uu,vv] = meshgrid(us,vs); %detector
    
    %Create weight according to each detector element
    w = (geo.DSD)./sqrt((geo.DSD)^2+uu.^2 + vv.^2);
    
    %Multiply the weights with projection data
    proj(:,:,ii) = proj(:,:,ii).*w';
end

%% filter
proj_filt = filtering(proj,geo,alpha); % Not sure if offsets are good in here
%RMFIELD Remove fields from a structure array.
geo=rmfield(geo,'filter');
%% backproject
res=Atb((proj_filt),geo,alpha); % Weighting is inside


if nargout>1
     error=proj-Ax(res,geo,alpha);
     errorL2=norm(error(:));
end

end