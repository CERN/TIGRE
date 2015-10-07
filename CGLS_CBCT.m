function x= CGLS_CBCT(proj,geo,angles,niter)

x=zeros(geo.nVoxel');
d=proj;
r0=Atb(proj,geo,angles);
t=Ax(r0,geo,angles);

p=r0;


for ii=1:niter
    alpha=norm(r0(:),2)^2/norm(t(:),2)^2
    x=x+alpha*p;
    d=d-alpha*t;
    r1=Atb(d,geo,angles);
    beta=norm(r1(:),2)^2/norm(r0(:),2)^2;
    p=r1+beta*p;
    t=Ax(p,geo,angles);
    
    r0=r1;
    
end




end