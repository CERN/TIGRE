function [x,errorL2]= CGLS_CBCT(proj,geo,angles,niter)

W=Ax(ones(geo.nVoxel'),geo,angles);  %
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;

% Back-Projection weigth, V
[x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
    -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
A = permute(angles, [1 3 2]);
V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
V=sum(V,3);
clear A x y dx dz;



reg=zeros(geo.nVoxel');

x=zeros(geo.nVoxel');
d=proj;
r0=bsxfun(@times,1./V,Atb(proj,geo,angles));
% r0=r0+0.01*reg;% Tikh Regularization;

t=W.*Ax(r0,geo,angles);
% reg=0.01*r0;% Tikh Regularization;



p=r0;

errorL2=zeros(niter,1);
for ii=1:niter
    
    alpha=norm(r0(:),2)^2/norm(t(:),2)^2;
    x=x+alpha*p;
    d=d-alpha*t;
    r1=bsxfun(@times,1./V,Atb(d,geo,angles));
%     r1=r1+0.01*reg;% Tikh Regularization;
    
    % Save error
    errorL2(ii)=norm(r1(:),2)^2;
    
    
    beta=norm(r1(:),2)^2/norm(r0(:),2)^2;
    p=r1+beta*p;
    t=W.*Ax(p,geo,angles);
%     reg=0.01*p; % Tikh Regularization;
    
    r0=r1;
    if mod(ii,50)==0
        disp(['Iteration: ',num2str(ii)]);
    end
    
end




end