
function [u] = TV_SB_3D_Time(f,geo, angles, N,mu, lambda, gamma, alpha, nInner, nBreg)


% Normalize data
% normFactor  = getNormalizationFactor(f,f);
% f           = normFactor*f;


if 1==2   %%%%%%%%%%%%% MANUCH
    % Normalize Jacobian such that its Hessian diagonal is equal to 1
    normFactorJ = 1/sqrt(max(diag(J'*J)));
    J           = J*normFactorJ;
    
    % Scale the forward and adjoint operations so doent depend on the size
    scale       = 1/max(abs(J'*f(:,1)));
    
    % Define forward and adjoint operators on each volume
    %A       = @(x)(((J*x(:)))/scale);
    %AT          = @(x)(reshape((J'*x)*scale,N(1:3)));
end
% Krylov convergence criterion: decrease to improve precision for solving
% the linear system, increase to go faster
tolKrylov   = 1e-2; % 1e-4

% Reserve memory for the auxillary variables
rows        = N(1);
cols        = N(2);
height      = N(3);
% time        = N(4);
f0          = f;
u           = zeros(N);
x           = zeros(N);
y           = zeros(N);
z           = zeros(N);
t           = zeros(N);
bx          = zeros(N);
by          = zeros(N);
bz          = zeros(N);
bt          = zeros(N);

murf=zeros([geo.nVoxel' size(f,4)],'single');
for it = 1:size(f,4)
    murf(:,:,:,it)=Atb(f(:,:,:,it),geo,angles,'matched');
end

%  Do the reconstruction
for outer = 1:nBreg
    for inner = 1:nInner
        % update u
        rhs     = murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by)+lambda*Dzt(z-bz) + ...
            + lambda*Dtt(t-bt);
        %u=reshape(krylov(rhs),N);
        
        u       = reshape(krylov(single(rhs(:))),N);   %%%%%%%%%%%%% MANUCH
        %      for i=1:size(f,4)
        %
        %          rhs1=rhs(:,:,:,time);
        %          im=CGLS(rhs1,geo,angles,3);
        %          u(:,:,:,time)=rhs;
        %      end
        
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
        %       fForw           = A(u(:,:,:,it));
        fForw=single(Ax(single(u(:,:,:,it)),geo,angles,'ray-voxel'));
        f(:,:,:,it)         = f(:,:,:,it) + f0(:,:,:,it) - fForw;
        %       murf(:,:,:,it)  = mu*AT(f(:,it));
        murf(:,:,:,it)=mu*single(Atb(f(:,:,:,it),geo,angles,'matched'));
    end
    
    
end


    function [xs,ys] = shrink2(x,y,lambda)
        s = sqrt(x.*conj(x)+y.*conj(y));
        ss = s-lambda;
        ss = ss.*(ss>0);
        
        s = s+(s<lambda);
        ss = ss./s;
        
        xs = ss.*x;
        ys = ss.*y;
    end

    function xs = shrink1(x,lambda)
        s = abs(x);
        xs = sign(x).*max(s-lambda,0);
    end

% =====================================================================
% Krylov solver subroutine
% X = GMRES(A,B,RESTART,TOL,MAXIT,M)
% bicgstab(A,b,tol,maxit)
    function dx = krylov(r)
        %            dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        %dx = gmres ((@jtjx), (r), []);
        [dx] = bicgstab(@jtjx, r, tolKrylov, 100);
        
    end

% =====================================================================
% Callback function for matrix-vector product (called by krylov)
    function b = jtjx(sol)
        
        solMat  = reshape(sol,N);
        
        % Laplacian part
        bTV     =single( lambda*(Dxt(Dx(solMat))+Dyt(Dy(solMat))+Dzt(Dz(solMat))+...
            Dtt(Dt(solMat))));
        
        
        % Jacobian part
        for it = 1:size(f,4)
            
            tmp             = single(solMat(:,:,:,it));
            %             bJac(:,:,:,it)  = mu*AT(A(tmp(:)));
            AAA=single(Ax(tmp,geo,angles,'ray-voxel'));
            bJac(:,:,:,it)=mu* single(Atb(AAA,geo,angles,'matched'));
            
        end
        
        % Stability term
        bG      = single(gamma*sol);
        
        b     = bTV(:) + bJac(:) + bG(:);
        
        
    end
% =====================================================================
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

