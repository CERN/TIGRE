function [u] = extended_split_bregman_TV(f,geo, angles,mu, lambda, gamma, alpha, nInner, nBreg,tolKrylov,max_iter)
% 4D split bregman TV reconstruction. From the following article
%
% Montesinos, P., Abascal, J.F.P., Cussó, L., Vaquero, J.J. and Desco, M. (2014),
% Application of the compressed sensing technique to self‐gated cardiac cine sequences in 
% small animals. Magn. Reson. Med., 
% 72: 369-380. https://doi.org/10.1002/mrm.24936

% Krylov convergence criterion: decrease to improve precision for solving
% the linear system, increase to go faster
if nargin<=9
    tolKrylov   = 1e-4;
    max_iter=100;
end
norm_factor=1./norm(f(:)/prod(geo.nDetector));
f=f.*norm_factor;

size4 = [geo.nVoxel',size(f,4)];
verbose=1;
% Reserve memory for the auxillary variables

f0          = f;
u           = zeros(size4,'single');
x           = zeros(size4,'single');
y           = zeros(size4,'single');
z           = zeros(size4,'single');
t           = zeros(size4,'single');
bx          = zeros(size4,'single');
by          = zeros(size4,'single');
bz          = zeros(size4,'single');
bt          = zeros(size4,'single');


murf=zeros([geo.nVoxel' size(f,4)],'single');
for it = 1:size(f,4)
    murf(:,:,:,it)=mu*Atb(f(:,:,:,it),geo,angles,'matched');
end

%  Do the reconstruction
for outer = 1:nBreg
    if (outer==1 && verbose==1);tic;end
    for inner = 1:nInner
        % update u
        rhs     = murf+lambda*(Dxt(x-bx)+Dyt(y-by)+Dzt(z-bz)+ Dtt(t-bt));
        u       = krylov(rhs,tolKrylov,max_iter,size4,lambda,mu,gamma,geo,angles);
        
        dx      = Dx(u);
        dy      = Dy(u);
        dz      = Dz(u);
        dt      = Dt(u);
        
        % update x and y and z
        [x,y]   = shrink2(dx+bx,dy+by,alpha(1)/lambda);
        z       = shrink1(dz+bz,alpha(2)/lambda);
        t       = shrink1(dt+bt,alpha(3)/lambda);
        
        % update bregman parameters
        bx          = bx+dx-x;
        by          = by+dy-y;
        bz          = bz+dz-z;
        bt          = bt+dt-t;
    end   % inner loop
    
    for it = 1:size(f,4)
        f(:,:,:,it)         = f(:,:,:,it) + f0(:,:,:,it) - Ax(single(u(:,:,:,it)),geo,angles,'ray-voxel');
        murf(:,:,:,it)=mu*single(Atb(f(:,:,:,it),geo,angles,'matched'));
    end
    if (outer==1 && verbose==1)
        expected_time=toc*nBreg;
        disp('extended split bregman TV');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
end
u=u./norm_factor;
end
%% Krylov solver subroutine
function dx = krylov(r,tolKrylov,max_iter,N,lambda,mu,gamma,geo,angles)

% alternatively.....
% dx = gmres (@(sol)jtjx(sol,N,lambda,mu,gamma,geo,angles), r(:), 30, tolKrylov, max_iter);
dx = bicgstab(@(sol)jtjx(sol,N,lambda,mu,gamma,geo,angles), r(:), tolKrylov, max_iter);
dx=reshape(dx,size(r));
end

% Callback function for matrix-vector product (called by krylov)
function b = jtjx(sol,N,lambda,mu,gamma,geo,angles)

sol  = reshape(single(sol),N);
% Laplacian part
bTV     = lambda*(Dxt(Dx(sol))+Dyt(Dy(sol))+Dzt(Dz(sol))+Dtt(Dt(sol)));
% Jacobian part
bJac=zeros(N,'single');
for it = 1:N(4)
    bJac(:,:,:,it)=mu* Atb(Ax(sol(:,:,:,it),geo,angles,'ray-voxel'),geo,angles,'matched');
end
% Stability term
b   = bTV(:) + bJac(:) + gamma*sol(:);


end
%% Shrinkage operators.
function [dx,dy] = shrink2(x,y,lambda)
s = sqrt(x.^2+y.^2);
ss = max(s-lambda,0)./(s+0.00001);

dx = ss.*x;
dy = ss.*y;
end
function xs = shrink1(x,lambda)
s = abs(x);
xs = sign(x).*max(s-lambda,0);
end
%% Numerical difference functions
% Backward differences in X
function d = Dy(u)
d=zeros(size(u));
d(:,2:end,:,:)=diff(u,1,2);
d(:,1,:,:)=u(:,1,:,:)-u(:,end,:,:);
end
% Transpose of the Backward differences in X
function d = Dyt(u)
d=zeros(size(u));
d(:,1:end-1,:,:)=-diff(u,1,2);
d(:,end,:,:)=(u(:,end,:,:)-u(:,1,:,:));
end
% Backward differences in Y
function d = Dx(u)
d=zeros(size(u));
d(2:end,:,:,:)=diff(u,1,1);
d(1,:,:,:)=u(1,:,:,:)-u(end,:,:,:);
end
% Transpose of the Backward differences in Y
function d = Dxt(u)
d=zeros(size(u));
d(1:end-1,:,:,:)=-diff(u,1,1);
d(end,:,:,:)=u(end,:,:,:)-u(1,:,:,:);
end
% Backward differences in Z
function d = Dz(u)
d=zeros(size(u));
d(:,:,2:end,:)=diff(u,1,3);
d(:,:,1,:)=u(:,:,1,:)-u(:,:,end,:);
end
% Transpose of the Backward differences in Z
function d = Dzt(u)
d=zeros(size(u));
d(:,:,1:end-1,:)=-diff(u,1,3);
d(:,:,end,:)=u(:,:,end,:)-u(:,:,1,:);
end
% Backward differences in T
function d = Dt(u)
d=zeros(size(u));
d(:,:,:,2:end)=diff(u,1,4);
d(:,:,:,1)=u(:,:,:,1)-u(:,:,:,end);
end
% Transpose of the Backward differences in T
function d = Dtt(u)
d=zeros(size(u));
d(:,:,:,1:end-1)=-diff(u,1,4);
d(:,:,:,end)=u(:,:,:,end)-u(:,:,:,1);
end

