function [x,errorL2]= CGLS_CBCT(proj,geo,angles,niter)



% //doi: 10.1088/0031-9155/56/13/004
x=zeros(geo.nVoxel'); % could be someting else

r=proj-Ax(x,geo,angles,'Krylov');
p=Atb(r,geo,angles,'Krylov');
gamma=norm(p(:),2)^2;


errorL2=zeros(niter,1);
for ii=1:niter
    
    
    q=Ax(p,geo,angles,'Krylov');
    alpha=gamma/norm(q(:),2)^2;
    x=x+alpha*p;
    r=r-alpha*q;
    
    s=Atb(r,geo,angles,'Krylov');
    gamma1=norm(s(:),2)^2;
    beta=gamma1/gamma;
    p=s+beta*p;
   
    if nargout>1
        aux=proj-Ax(x,geo,angles,'Krylov');
        errorL2(ii)=norm(aux(:),2);
    end
    if mod(ii,50)==0
        disp(['Iteration: ',num2str(ii)]);
    end
    
end




end