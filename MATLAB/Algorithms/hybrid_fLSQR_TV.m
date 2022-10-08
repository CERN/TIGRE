function [x,errorL2,lambda_vec, qualMeasOut]= hybrid_fLSQR_TV(proj,geo,angles,niter,varargin)

% LSQR solves the CBCT problem using LSQR. 
% This is mathematically equivalent to CGLS.
% 
%  LSQR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
% 
%  LSQR(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
% 
% 
%  'Init'    Describes diferent initialization techniques.
%             * 'none'     : Initializes the image to zeros (default)
%             * 'FDK'      : intializes image to FDK reconstrucition
%             * 'multigrid': Initializes image by solving the problem in
%                            small scale and increasing it when relative
%                            convergence is reached.
%             * 'image'    : Initialization using a user specified
%                            image. Not recomended unless you really
%                            know what you are doing.
%  'InitImg'    an image for the 'image' initialization. Avoid.
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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:           max(K)   https://github.com/CERN/TIGRE/
% Coded by:           Malena Sabate Landman, Ander Biguri
%--------------------------------------------------------------------------

%%

[verbose,x0,QualMeasOpts,gpuids, lambda, NoiseLevel]=parse_inputs(proj,geo,angles,varargin);

%%% PARAMETER CHOICE HIERARCHY: given lambda, DP, GCV

if isnan(lambda)
    if isnan(NoiseLevel)
        RegParam = 'gcv';
        % Malena: if this is not good, we can use alternative formulations
    else
        RegParam = 'discrepit'; 
        % Malena: if this is slow, we can use an adaptive method 
    end
else
    RegParam = 'given_lambda';
end

% msl: no idea of what this is. Should I check?
measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter);

% Paige and Saunders //doi.org/10.1145/355984.355989

% Initialise matrices
U = single(zeros(prod(size(proj)), niter+1));
V = single(zeros(prod(geo.nVoxel), niter)); % Malena: Check if prod(geo.nVoxel) is correct, I want size of object
Z = single(zeros(prod(geo.nVoxel), niter)); % Flexible basis
M = zeros(niter+1,niter); % Projected matrix 1
T = zeros(niter+1); % Projected matrix 2
proj_rhs = zeros(niter+1,1); % Projected right hand side
lambda_vec = zeros(niter,1);

z = zeros(geo.nVoxel');

% Enumeration as given in the paper for 'Algorithm LSQR'
% (1) Initialize 
% norm_proj = norm (proj(:),2);
u = proj-Ax(x0,geo,angles,'Siddon','gpuids',gpuids);

% Build xA0 
null_D = single((prod(geo.nVoxel)^(-1/2))*ones(geo.nVoxel'));

k_aux = Ax(null_D,geo,angles,'Siddon','gpuids',gpuids);
xA0 = (prod(geo.nVoxel)^(-1/2))*norm(k_aux(:),2)^(-2)* (k_aux(:)'*u(:)); % scalar
xA0 = single(ones(geo.nVoxel')) * xA0;

% Other way of computing the same, what is better?
% vecAN = Ax(null_D,geo,angles,'Siddon','gpuids',gpuids);
% [Q0, R0] = qr(vecAN(:), 0);
% x0 = Q0'*u(:); x0 = R0\x0; x0 = null_D(:)*x0;
% Ax0 = A_times_vec(A, x0)

k_aux = (prod(geo.nVoxel)^(-1/2))*norm(k_aux(:),2)^(-2)*Atb(k_aux,geo,angles,'matched','gpuids',gpuids);

u = u-Ax(xA0,geo,angles,'Siddon','gpuids',gpuids);

normr = norm(u(:),2);
u = u/normr;
U(:,1) = u(:);

if max(max(max(x0))) == 0
    W = ones(size(x0));
else
    W = build_weights (x0);
end

proj_rhs(1) = normr;
errorL2 = zeros(1,niter);


% (2) Start iterations 
for ii=1:niter

    if (ii==1 && verbose);tic;end

    % Update V, Z, and projected matrix T
    v = Atb(u,geo,angles,'matched','gpuids',gpuids);
    for jj = 1:ii-1
        T(jj,ii) = V(:,jj)'*v(:);
        v(:) = v(:) - T(jj,ii) *V(:,jj);
    end
    T(ii,ii) = norm(v(:),2);
    v = v/T(ii,ii);
    % Choose if fixing a tolerance or a number of iterations of both
    % This should maybe be done in single precisions...
    z(:) = mvpE(k_aux, v(:), 'transp');
    aux_z = lsqr(@(x,tflag) Ltx(W, x, tflag), z(:), [], 50);
    z(:) = lsqr(@(x,tflag) Lx(W, x, tflag), aux_z, [], 50);
    z(:) = mvpE(k_aux, z(:), 'notransp');
%     z(:) = lsqr(@(x,tflag)mvpEt(k_aux, x, tflag), v(:), [], 25);
%     z(:) = lsqr(@(x,tflag)mvpE(k_aux, x, tflag), z(:), [], 25);
    z = single(z);

    V(:,ii) = v(:);
    Z(:,ii) = z(:);

    % Update U and projected matrix M
    u = Ax(z,geo,angles,'Siddon','gpuids',gpuids);
    for jj = 1:ii
        M(jj,ii) = U(:,jj)'*u(:);
        u(:) = u(:) - M(jj,ii) *U(:,jj);
    end
    M(ii+1,ii) = norm(u(:),2);
    u = u / M(ii+1,ii);
    U(:,ii+1) = u(:);

    % Malena. Proposed update: we should check algorithms breaks; 
    % 'if abs(alpha) <= eps || abs(beta) <= eps' - end and save

    %%% Solve the regularised projected problem 
    % (using the SVD of the small projected matrix)
    Mk = M(1:ii+1,1:ii);

    % Prepare the projected regularization term
    WZ = zeros(3*prod(geo.nVoxel),ii);
    for jj=1:ii
        % This can be done more efficiently... 
        % DZ can be saved and updates at each iteration
        out = Lx (W, Z(:,jj), 'notransp');
        WZ(:,jj) = out; %[out{1}(:);out{2}(:);out{3}(:)];
    end
    [~, ZRk] = qr(WZ(:,1:ii), 0);
    ZRksq = ZRk(1:ii,1:ii);
    rhsk = proj_rhs(1:ii+1);
    
    if strcmp(RegParam,'discrepit')
        eta = 1.01;
        if discrepancy_Tik(0, Mk, ZRksq, rhsk, eta*NoiseLevel) > 0
            lambda = 0;
        else
            lambda = fzero(@(l)discrepancy_Tik(l, Mk, ZRksq, rhsk, eta*NoiseLevel), [0, 1e10]);
        end
        lambda_vec(ii) = lambda; % We should output this, maybe?
    elseif strcmp(RegParam,'gcv')
        [Uk,  ~, ~, Ck, Sk] = gsvd(Mk, ZRksq);
        rhskhat = Uk'*rhsk;
        if ii==1
            gammak = Ck(1)/Sk(1);
        else
            gammak = sqrt(diag(Ck'*Ck)./diag(Sk'*Sk));
        end
        lambda = fminbnd(@(l)gcv(l, rhskhat, gammak),  0, double(gammak(ii)));
        lambda_vec(ii) = lambda; % We should output this, maybe?

    elseif strcmp(RegParam,'given_lambda')
        lambda_vec(ii) = lambda;
    end

    MZk = [Mk; lambda*ZRksq];
    rhsZk = [rhsk; zeros(ii,1)];
    y = MZk\rhsZk;

    errorL2(ii)=norm(rhsk - Mk*y);

    d = Z(:,1:ii)*y;
    x = x0 + reshape(d,size(x0)) + xA0;

    W = build_weights (x);

    % Test for convergence. 
    % msl: I still need to implement this. 
    % msl: There are suggestions on the original paper. Let's talk about it!
    
    if measurequality 
        qualMeasOut(:,ii)=Measure_Quality(x0,x,QualMeasOpts);
    end

    if (ii==1 && verbose)
        expected_time=toc*niter;   
        disp('LSQR');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);   
        disp('');
     end
end

end

%%% Regularization parameter choices

function out = discrepancy_Tik(lambda, A, ZRksq, b, nnoise)
n = size(A,2);
if n == 1
    out = 0;
else
    xl = [A; lambda*ZRksq]\[b; zeros(size(ZRksq,1),1)];
    out = (norm(A*xl - b)/norm(b))^2 - nnoise^2;
end
end 

function out = gcv(lambda, bhat, s)
% GCV for the projected problem - no weights
% If Bk is the projected matrix and Bk=Uk*Sk*Vk^T
% lambda is the regularisation parameter
% bhat is Uk'*bk 
% s=diag(Sk) 

m = length(bhat);
n = length(s);

t0 = sum(abs(bhat(n+1:m)).^2);

s2 = abs(s) .^ 2;
lambda2 = lambda^2;

t1 = lambda2 ./ (s2 + lambda2);
t2 = abs(bhat(1:n) .* t1) .^2;

out = (sum(t2) + t0) / ((sum(t1)+m-n)^2);

end 

%%%%%%

function W = build_weights ( x)
    % Directional discrete derivatives 
    % Reflective boundary conditions 
    Dxx=zeros(size(x),'single');
    Dyx=zeros(size(x),'single');
    Dzx=zeros(size(x),'single');
    
    Dxx(1:end-1,:,:)=x(1:end-1,:,:)-x(2:end,:,:);
    Dyx(:,1:end-1,:)=x(:,1:end-1,:)-x(:,2:end,:);
    Dzx(:,:,1:end-1)=x(:,:,1:end-1)-x(:,:,2:end);

    W = (Dxx.^2+Dyx.^2+Dzx.^2+1e-6).^(-1/4); % Fix this...

end

function out = Lx (W, x, transp_flag)
if strcmp(transp_flag,'transp')
    x = reshape(x,[size(W),3]);
    Wx_1 = W .* x(:,:,:,1);
    Wx_2 = W .* x(:,:,:,2);
    Wx_3 = W .* x(:,:,:,3);

    % Left here, but this is how to make a transpose
    DxtWx_1=Wx_1;
    DytWx_2=Wx_2;
    DztWx_3=Wx_3;
    
    DxtWx_1(2:end-1,:,:)=Wx_1(2:end-1,:,:)-Wx_1(1:end-2,:,:);
    DxtWx_1(end,:,:)=-Wx_1(end-1,:,:);
    
    DytWx_2(:,2:end-1,:)=Wx_2(:,2:end-1,:)-Wx_2(:,1:end-2,:);
    DytWx_2(:,end,:)=-Wx_2(:,end-1,:);
    
    DztWx_3(:,:,2:end-1)=Wx_3(:,:,2:end-1)-Wx_3(:,:,1:end-2);
    DztWx_3(:,:,end)=-Wx_3(:,:,end-1);

    out = DxtWx_1 + DytWx_2 + DztWx_3;
    out = out(:);

elseif strcmp(transp_flag,'notransp')
    
    x = reshape(x,size(W));

    % Directional discrete derivatives
    % Reflective boundary conditions 
	Dxx=zeros(size(x),'single');
    Dyx=zeros(size(x),'single');
    Dzx=zeros(size(x),'single');
    
    Dxx(1:end-1,:,:)=x(1:end-1,:,:)-x(2:end,:,:);
    Dyx(:,1:end-1,:)=x(:,1:end-1,:)-x(:,2:end,:);
    Dzx(:,:,1:end-1)=x(:,:,1:end-1)-x(:,:,2:end);
    % Build weights - is it better to find the right rotation and add
    % tensors?
    out = cat(4,W .* Dxx, W .* Dyx, W .* Dzx);
    out = out(:);
end
end

function out = Ltx (W, x, transp_flag)
if strcmp(transp_flag,'notransp')
    x = reshape(x,[size(W),3]);
    Wx_1 = W .* x(:,:,:,1);
    Wx_2 = W .* x(:,:,:,2);
    Wx_3 = W .* x(:,:,:,3);

    % Left here, but this is how to make a transpose
    DxtWx_1=Wx_1;
    DytWx_2=Wx_2;
    DztWx_3=Wx_3;
    
    DxtWx_1(2:end-1,:,:)=Wx_1(2:end-1,:,:)-Wx_1(1:end-2,:,:);
    DxtWx_1(end,:,:)=-Wx_1(end-1,:,:);
    
    DytWx_2(:,2:end-1,:)=Wx_2(:,2:end-1,:)-Wx_2(:,1:end-2,:);
    DytWx_2(:,end,:)=-Wx_2(:,end-1,:);
    
    DztWx_3(:,:,2:end-1)=Wx_3(:,:,2:end-1)-Wx_3(:,:,1:end-2);
    DztWx_3(:,:,end)=-Wx_3(:,:,end-1);

    out = DxtWx_1 + DytWx_2 + DztWx_3;
    out = out(:);

elseif strcmp(transp_flag,'transp')
    
    x = reshape(x,size(W));

    % Directional discrete derivatives
    % Reflective boundary conditions 
	Dxx=zeros(size(x),'single');
    Dyx=zeros(size(x),'single');
    Dzx=zeros(size(x),'single');
    
    Dxx(1:end-1,:,:)=x(1:end-1,:,:)-x(2:end,:,:);
    Dyx(:,1:end-1,:)=x(:,1:end-1,:)-x(:,2:end,:);
    Dzx(:,:,1:end-1)=x(:,:,1:end-1)-x(:,:,2:end);
    % Build weights - is it better to find the right rotation and add
    % tensors?
    out = cat(4,W .* Dxx, W .* Dyx, W .* Dzx);
    out = out(:);
end
end

function out = mvpE(k_aux, x , transp_flag)
if strcmp(transp_flag,'transp')
    out = x(:) - k_aux(:)*(ones(size(x(:)))'*x(:));
elseif strcmp(transp_flag,'notransp')
    out = x(:) - ones(size(x(:)))*(k_aux(:)'*x(:));
end
end

% function out = mvpEt(k_aux, x , transp_flag)
% if strcmp(transp_flag,'notransp')
%     out = x(:) - k_aux(:)*(ones(size(x(:)))'*x(:));
% elseif strcmp(transp_flag,'transp')
%     out = x(:) - ones(size(x(:)))*(k_aux(:)'*x(:));
% end
% end


%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids, lambda, NoiseLevel]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','lambda','noiselevel'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:LSQR:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
       error('TIGRE:LSQR:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]); 
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
   if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end
        if isempty(ind)
            error('TIGRE:LSQR:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
        end
        val=argin{jj};
    end
    
    switch opt
        case 'init'
            x=[];
            if default || strcmp(val,'none')
                x=zeros(geo.nVoxel','single');
                continue;
            end
            if strcmp(val,'FDK')
                x=FDK(proj,geo,angles);
                continue;
            end
            if strcmp(val,'multigrid')
                x=init_multigrid(proj,geo,angles);
                continue;
            end
            if strcmp(val,'image')
                initwithimage=1;
                continue;
            end
            if isempty(x)
               error('TIGRE:LSQR:InvalidInput','Invalid Init option') 
            end
            % % % % % % % ERROR
        case 'initimg'
            if default
                continue;
            end
            if exist('initwithimage','var')
                if isequal(size(val),geo.nVoxel')
                    x=single(val);
                else
                    error('TIGRE:LSQR:InvalidInput','Invalid image for initialization');
                end
            end
        case 'lambda'
            if default
                lambda=NaN;
            else
                lambda=val;
            end
        case 'noiselevel'
            if default
                NoiseLevel=NaN;
            else
                NoiseLevel=val;
            end

        %  =========================================================================
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:LSQR:InvalidInput','Invalid quality measurement parameters');
                end
            end
         case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE:LSQR:Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
        otherwise 
            error('TIGRE:LSQR:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end