function [res,errorL2,rmFDK,corrFDK,mssFDK]=FDK_CBCT(proj,geo,alpha,p,filter,varargin)
%Input is data,geosize,angles

if nargin<5
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
res=Atb(proj_filt,geo,alpha); % Weighting is inside


rmFDK=RMSE(p,res);
corrFDK=CC(p,res);
mssFDK=MSSIM(p,res);

if nargout>1
     error=proj-Ax(res,geo,alpha);
     errorL2=norm(error(:));
end

end