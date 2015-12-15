function [x,errorL2]= CGLS_CBCT(proj,geo,angles,niter)



% //doi: 10.1088/0031-9155/56/13/004
x=zeros(geo.nVoxel'); % could be someting else

r=proj-Ax(x,geo,angles,'Krylov');
p=Atb(r,geo,angles,'Krylov');
gamma=norm(p(:),2)^2;


errorL2=zeros(niter,1);
for ii=1:niter
     if ii==1;tic;end
    
    q=Ax(p,geo,angles,'Krylov');
    alpha=gamma/norm(q(:),2)^2;
    x=x+alpha*p;
    r=r-alpha*q;
    
    s=Atb(r,geo,angles,'Krylov');
    gamma1=norm(s(:),2)^2;
    beta=gamma1/gamma;
    gamma=gamma1;
    p=s+beta*p;
   
    % Diverges and that is not cool. I dont know why. Paper says tehre is
    % less than 1% of error between the A and At, but could that bee too
    % much anyway?
    if nargout>1
        aux=proj-Ax(x,geo,angles,'Krylov');
        errorL2(ii)=norm(aux(:),2);
    end
    if ii>1 && errorL2(ii)>errorL2(ii-1)
       x=x-alpha*p;
       return; 
    end
     if ii==1;
        expected_time=toc*niter;   
        disp('CGLS');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);   
        disp('');
    end
end




end