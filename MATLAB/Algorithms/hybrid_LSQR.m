function [x,resL2, qualMeasOut,lambda_vec]= hybrid_LSQR(proj,geo,angles,niter,varargin)

% hybrid_LSQR solves the CBCT problem using LSQR.
%
%  hybrid_LSQR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%  hybrid_LSQR(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%  'lambda'  Value of parameter lambda, default autocomputed.
%  'Noiselevel' the expected nosie level, in %, replaces lambda.
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
%  'groundTruth'  an image as grounf truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
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
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Malena Sabate Landman, Ander Biguri
%--------------------------------------------------------------------------

%%

[verbose,x0,QualMeasOpts,gpuids, lambda, NoiseLevel,gt]=parse_inputs(proj,geo,angles,varargin);

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
measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    x_prev=gt;
    clear gt
end
if nargout<3 && measurequality
    warning("Image metrics requested but none catched as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),niter);
resL2 = zeros(1,niter);
% Paige and Saunders //doi.org/10.1145/355984.355989

% This should be one to avoid loosing orthogonality, but can be switched
% off for testing.
reorth = 1;

% Initialise matrices
U = zeros(numel(proj), niter+1,'single');
V = zeros(prod(geo.nVoxel), niter,'single'); % Malena: Check if prod(geo.nVoxel) is correct
B = zeros(niter+1,niter,'single'); % Projected matrix
proj_rhs = zeros(niter+1,1,'single'); % Projected right hand side
lambda_vec = zeros(niter,1);

% Enumeration as given in the paper for 'Algorithm LSQR'
% (1) Initialize
u = proj-Ax(x0,geo,angles,'Siddon','gpuids',gpuids);
normr = norm(u(:),2);
u = u/normr;
U(:,1) = u(:);

beta = normr;
proj_rhs(1) = normr;

% (2) Start iterations
for ii=1:niter
    if measurequality && ~strcmp(QualMeasOpts,'error_norm')
        if ii==1
            x_prev=x0;
        else
            x_prev =  x0 + reshape(V(:,1:ii)*y,size(x0)); % only store if necesary
        end
    end
    if (ii==1 && verbose);tic;end
    
    % Update V_ii
    v = Atb(u,geo,angles,'matched','gpuids',gpuids);
    
    if ii>1
        v(:) = v(:) - beta*V(:,ii-1);
    end
    if reorth % Maybe change to matrix operations!
        for jj = 1:ii-1
            v(:) = v(:) - (V(:,jj)'*v(:))*V(:,jj);
        end
    end
    alpha = norm(v(:),2); % msl: do we want to check if it is 0?
    v = v/alpha;
    V(:,ii) = v(:);
    
    % Update U_{ii+1}
    u = Ax(v,geo,angles,'Siddon','gpuids',gpuids) - alpha*u;
    if reorth % Maybe change to matrix operations!
        for jj = 1:ii-1
            u(:) = u(:) - (U(:,jj)'*u(:))*U(:,jj);
        end
    end
    beta = norm(u(:),2);
    u = u / beta;
    U(:,ii+1) = u(:);
    
    % Update projected matrix
    B(ii,ii) = alpha;
    B(ii+1,ii) = beta;
    % Malena. Proposed update: we should check algorithms breaks;
    % 'if abs(alpha) <= eps || abs(beta) <= eps' - end and save
    
    % Solve the projected problem
    % (using the SVD of the small projected matrix)
    Bk = B(1:ii+1,1:ii);
    [Uk, Sk, Vk] = svd(Bk);
    if ii==1
        Sk = Sk(1,1);
    else
        Sk = diag(Sk);
    end
    rhsk = proj_rhs(1:ii+1);
    rhskhat = Uk'*rhsk;
    lsqr_res = abs(rhskhat(ii+1))/normr;
    
    
    if strcmp(RegParam,'discrepit')
        eta = 1;
        if lsqr_res > eta*NoiseLevel
            lambda = 0;
        else
            lambda = fzero(@(l)discrepancy(l, Bk, rhsk, eta*NoiseLevel), [0, 1e10]);
        end
        lambda_vec(ii) = lambda; % We should output this, maybe?
    elseif strcmp(RegParam,'gcv')
        lambda = fminbnd(@(l)gcv(l, rhskhat, Sk),  0, double(Sk(1)));
        % Adapt from IRtools
        lambda_vec(ii) = lambda; % We should output this, maybe?
    elseif strcmp(RegParam,'given_lambda')
        lambda_vec(ii) = lambda;
    end
    
    Dk = Sk.^2 + lambda^2;
    rhskhat = Sk .* rhskhat(1:ii);
    yhat = rhskhat(1:ii)./Dk;
    y = Vk * yhat;
    
%     resL2(ii)=norm(rhsk - Bk*y); % residual norm
    x = x0 + reshape(V(:,1:ii)*y,size(x0));
    if measurequality
        qualMeasOut(:,ii)=Measure_Quality(x_prev,x,QualMeasOpts);
    end
    aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
    resL2(ii)=im3Dnorm(aux,'L2');
    if ii>1 && resL2(ii)>resL2(ii-1)
        disp(['Algorithm stoped in iteration ', num2str(ii),' due to loss of ortogonality.'])
        return;
    end
    
    if (ii==1 && verbose)
        expected_time=toc*niter;
        disp('hybrid LSQR');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
end
x = x0 + reshape(V(:,1:ii)*y,size(x0));


end

%%% Regularization parameter choices

function out = discrepancy(l, A, b, nnoise)
n = size(A,2);
if n == 1
    out = 0;
else
    xl = [A; l*eye(n)]\[b; zeros(n,1)];
    out = (norm(A*xl -b)/norm(b))^2 - nnoise^2;
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


%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids,lambda,NoiseLevel,gt]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','lambda','noiselevel','groundtruth'};
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
        case 'groundtruth'
            if default
                gt=nan;
            else
                gt=val;
            end
        otherwise
            error('TIGRE:LSQR:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end