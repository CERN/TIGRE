function [tvg] =weighted_gradientTVnorm(img,delta)

tvg=zeros(size(img));

%Compute the local image-intensity

Gx=diff(img,1,1);
Gy=diff(img,1,2);
Gz=diff(img,1,3);


Gx=cat(1,zeros(size(Gx(1,:,:))),Gx);
Gy=cat(2,zeros(size(Gy(:,1,:))),Gy);
Gz=cat(3,zeros(size(Gz(:,:,1))),Gz);



%--------------------------------------------------------------------------
%Weighting equation#1 (An exponential function of the component of the local image-intensity gradient vector)


wx=exp(-(Gx./delta).^2);
wy=exp(-(Gy./delta).^2);
wz=exp(-(Gz./delta).^2);
%--------------------------------------------------------------------------

nrm=weighted_safenorm(Gx,Gy,Gz,wx,wy,wz);

tvg(1:end,1:end,1:end)= tvg(1:end,1:end,1:end)+((wx(1:end,1:end,1:end).*Gx(1:end,1:end,1:end))+(wy(1:end,1:end,1:end).*Gy(1:end,1:end,1:end))+(wz(1:end,1:end,1:end).*Gz(1:end,1:end,1:end)))./nrm(1:end,1:end,1:end);
tvg(2:end-1,:,:)=tvg(2:end-1,:,:)-(wx([2:end-1]+1,:,:).*Gx([2:end-1]+1,:,:))./nrm([2:end-1]+1,:,:);
tvg(:,2:end-1,:)=tvg(:,2:end-1,:)-(wy(:,[2:end-1]+1,:).*Gy(:,[2:end-1]+1,:))./nrm(:,[2:end-1]+1,:);
tvg(:,:,2:end-1)=tvg(:,:,2:end-1)-(wz(:,:,[2:end-1]+1).*Gz(:,:,[2:end-1]+1))./nrm(:,:,[2:end-1]+1);




end

function nrm=weighted_safenorm(Gx,Gy,Gz,wx,wy,wz)

ksi=0.00000001;
nrm=sqrt(ksi+(wx.*Gx.^2)+(wy.*Gy.^2)+(wz.*Gz.^2));
nrm(nrm==0)=1;

end

