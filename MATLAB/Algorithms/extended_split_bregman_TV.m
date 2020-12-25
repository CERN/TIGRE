function [u] = extended_split_bregman_TV(f,geo, angles, N,mu, lambda, gamma, alpha, nInner, nBreg)
% 4D split bregman TV reconstruction. From the following article
%
% Montesinos, P., Abascal, J.F.P., Cussó, L., Vaquero, J.J. and Desco, M. (2014),
% Application of the compressed sensing technique to self‐gated cardiac cine sequences in 
% small animals. Magn. Reson. Med., 
% 72: 369-380. https://doi.org/10.1002/mrm.24936

% Krylov convergence criterion: decrease to improve precision for solving
% the linear system, increase to go faster
tolKrylov   = 1e-4; % 1e-4

% Reserve memory for the auxillary variables

f0          = f;
u           = zeros(N,'single');
x           = zeros(N,'single');
y           = zeros(N,'single');
z           = zeros(N,'single');
t           = zeros(N,'single');
bx          = zeros(N,'single');
by          = zeros(N,'single');
bz          = zeros(N,'single');
bt          = zeros(N,'single');
max_iter=100;

murf=zeros([geo.nVoxel' size(f,4)],'single');
for it = 1:size(f,4)
    murf(:,:,:,it)=Atb(f(:,:,:,it),geo,angles,'matched');
end

%  Do the reconstruction
for outer = 1:nBreg
    for inner = 1:nInner
        % update u
        rhs     = murf+lambda*(Dxt(x-bx)+Dyt(y-by)+Dzt(z-bz)+ Dtt(t-bt));
        u       = krylov(rhs,tolKrylov,max_iter,N,lambda,mu,gamma,geo,angles);
        
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
    
end

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

solMat  = reshape(single(sol),N);
% Laplacian part
bTV     = lambda*(Dxt(Dx(solMat))+Dyt(Dy(solMat))+Dzt(Dz(solMat))+Dtt(Dt(solMat)));
% Jacobian part
bJac=zeros(N,'single');
for it = 1:N(4)
    bJac(:,:,:,it)=mu* Atb(Ax(solMat(:,:,:,it),geo,angles,'ray-voxel'),geo,angles,'matched');
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
function d = Dx(u)
% TODO: I think  this is doing the wrong axes (i.e. this is Dy)
d=diff(u,1,2);
d=cat(2,u(:,1,:,:)-u(:,end,:,:),d);
end
% Transpose of the Backward differences in X
function d = Dxt(u)
% TODO: I think  this is doing the wrong axes (i.e. this is Dy)
d=-diff(u,1,2);
d=cat(2,d,(u(:,end,:,:)-u(:,1,:,:)));
end
% Backward differences in Y
function d = Dy(u)
% TODO: I think  this is doing the wrong axes (i.e. this is Dx)
d=diff(u,1,1);
d=cat(1,u(1,:,:,:)-u(end,:,:,:),d);
end
% Transpose of the Backward differences in Y
function d = Dyt(u)
% TODO: I think  this is doing the wrong axes (i.e. this is Dx)
d=-diff(u,1,1);
d=cat(1,d,(u(end,:,:,:)-u(1,:,:,:)));
end
% Backward differences in Z
function d = Dz(u)
d=diff(u,1,3);
d=cat(3,u(:,:,1,:)-u(:,:,end,:),d);
end
% Transpose of the Backward differences in Z
function d = Dzt(u)
d=-diff(u,1,3);
d=cat(3,d,(u(:,:,end,:)-u(:,:,1,:)));
end
% Backward differences in T
function d = Dt(u)
d=diff(u,1,4);
d=cat(4,u(:,:,:,1)-u(:,:,:,end),d);
end
% Transpose of the Backward differences in T
function d = Dtt(u)
d=-diff(u,1,4);
d=cat(4,d,(u(:,:,:,end)-u(:,:,:,1)));
end

