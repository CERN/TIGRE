function [x_outer,resL2,qualMeasOut]= IRN_PICCS_CGLS(proj,geo,angles,niter,prior,varargin)
% IRN_PICCS_CGLS solves the CBCT problem with prior image constrained
% compressed sensing (PICCS), by iteratively reweighted norms and CGLS
%
%  IRN_PICCS_CGLS(PROJ,GEO,ANGLES,NITER,PRIOR) solves
%
%     min_x ||Ax - b||_2^2 + lambda*TV(x) + lambda_prior*TV(x - PRIOR)
%
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry described in GEO, using NITER CGLS iterations in total.
%   PRIOR is an earlier image of the same object, of size geo.nVoxel', for
%   example from a previous, fully sampled scan. Each TV term is replaced by
%   a weighted least-squares term, with weights rebuilt from the current
%   image at each of 'niter_outer' outer iterations, and each weighted
%   problem is solved by CGLS, as in IRN_TV_CGLS.
%
%  IRN_PICCS_CGLS(PROJ,GEO,ANGLES,NITER,PRIOR,OPT,VAL,...) uses options and
%   values for solving. The possible options in OPT are:
%
%  'lambda'       Weight of the TV(x) term, default 10.
%  'lambda_prior' Weight of the prior term, default 10.
%  'priorNorm'    The norm of the prior term:
%             * 'tv' : TV(x - PRIOR), PICCS (default)
%             * 'l2' : ||x - PRIOR||_2^2, as in PIPLE
%             * 'l1' : ||x - PRIOR||_1
%  'niter_outer'  Number of outer (reweighting) iterations, default 10.
%                 NITER is shared between them.
%  'outerRestart' What each outer iteration starts from:
%             * 'cold' : the initial image (default)
%             * 'warm' : the image of the previous outer iteration
%  'weights'      Weights for the TV(x) term in the first outer iteration,
%                 of size geo.nVoxel', instead of computing them from the
%                 initial image (for example from a previous solve).
%  'weightsPrior' Same, for the prior term ('tv' and 'l1' only).
%  'Init'    Describes different initialization techniques.
%             * 'none'     : Initializes the image to zeros (default)
%             * 'FDK'      : initializes image to FDK reconstruction
%             * 'multigrid': Initializes image by solving the problem in
%                            small scale and increasing it when relative
%                            convergence is reached.
%             * 'image'    : Initialization using a user specified
%                            image. Not recommended unless you really
%                            know what you are doing.
%  'InitImg'    an image for the 'image' initialization. Avoid.
%  'groundTruth'  an image as ground truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
%  'restart'  true or false. By default, when orthogonality is lost CGLS
%             restarts once from the current image; if that does not help
%             (or with false) the outer iteration ends and the weights are
%             rebuilt.
%  'checkResidual' true or false. By default the algorithm computes the
%             residual every iteration, which costs one extra forward
%             projection, and watches the residual of the whole regularised
%             system (data, TV and prior terms), which CGLS never lets grow:
%             if it grows, orthogonality was lost. false skips the check,
%             for speed; ||Ax-b|| is still computed if asked as an output.
%
% OUTPUTS:
%
%    [img]                       will output the reconstructed image
%    [img,resL2]                 will also output the L2 norm of the
%                                residual, ||Ax-b||, after each iteration
%    [img,resL2,qualMeas]        will also output the quality measurements
%                                asked by the input 'QualMeas'
%
%  F. Hastings, S. M. R. S. Islam, M. Sabate Landman, S. Hatamikia,
%  C.-B. Schonlieb, A. Biguri. Real-time CBCT reconstructions using Krylov
%  solvers in repeated scanning procedures. Phys. Med. Biol. 71 (2026)
%  105019. https://doi.org/10.1088/1361-6560/ae676e
%
%  See also IRN_TV_CGLS, PICCS
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
% Coded by:           Malena Sabate Landman, Ander Biguri, Fred Hastings,
%                     S M Ragib Shahriar Islam
%--------------------------------------------------------------------------

% An Iteratively Reweighted Norm Algorithm for Total Variation Regularization
% Paul Rodriguez and Brendt Wohlberg
% 10.1109/ACSSC.2006.354879

[verbose,x0,QualMeasOpts,gpuids,lambda,lambda_prior,priorNorm,niter_outer,gt,restart,checkResidual,outerRestart,W1_init,W2_init]=parse_inputs(proj,geo,angles,varargin);
% IRN_TV_CGLS calls this with an empty prior
usePrior=~isempty(prior) && lambda_prior>0;
if isempty(prior)
    name='IRN_TV_CGLS';
else
    name='IRN_PICCS_CGLS';
    if ~isequal(size(prior),geo.nVoxel')
        error('TIGRE:IRN_PICCS_CGLS:InvalidInput','The prior should be of size geo.nVoxel''');
    end
    prior=single(prior);
end
computeResidual=checkResidual || nargout>1;
x=x0;
x_outer=x0;

niter_break=round(niter/niter_outer);

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
iter=0;
sl=sqrt(lambda);
slp=sqrt(lambda_prior);
W2=[];
while iter<niter
    if (iter==0 && verbose);tic;end

    % Weights from the current image, N1 x N2 x N3 (or the ones given, the first time)
    if ~isempty(W1_init)
        W1=W1_init;
        W1_init=[];
    else
        W1 = build_weights(x);
    end
    if usePrior && ~strcmp(priorNorm,'l2')
        if ~isempty(W2_init)
            W2=W2_init;
            W2_init=[];
        elseif strcmp(priorNorm,'tv')
            W2 = build_weights(x-prior);
        else
            W2 = ((x-prior).^2+1e-6).^(-1/4); % ||W2.*v||^2 is a surrogate of ||v||_1
        end
    end

    if strcmp(outerRestart,'cold')
        x=x0;
    end

    % CGLS on the stacked system [A; sl*W1*D; slp*P] x = [b; 0; slp*P*prior].
    % Loss of orthogonality is detected on the residual of that whole system,
    % which CGLS keeps non-increasing in exact arithmetic. ||Ax-b|| alone is
    % allowed to rise when the regularisation is strong, so it is not used.
    restartCGLS=true;
    fullres_prev=inf;
    restarted_at=-1;
    iter_outer_start=iter;
    while iter<niter %this will be broken by a parameter
        if restartCGLS
            r_aux_1 = proj - Ax(x,geo,angles,'Siddon','gpuids',gpuids);
            r_aux_2 = -sl*Lx(W1, x);
            p = Atb(r_aux_1 ,geo,angles,'matched','gpuids',gpuids) + sl*Ltx(W1, r_aux_2);
            if usePrior
                r_aux_3 = -slp*Px(W2, x-prior, priorNorm);
                p = p + slp*Ptx(W2, r_aux_3, priorNorm);
            end
            gamma=norm(p(:),2)^2;
            restartCGLS=false;
        end
        iter=iter+1;
        if measurequality && ~any(strcmp(QualMeasOpts,'error_norm'))
            x_prev = x; % only store if necessary
        end
        q_aux_1 = Ax(p ,geo,angles,'Siddon','gpuids',gpuids);
        q_aux_2 = sl*Lx(W1, p);
        qq = norm(q_aux_1(:))^2 + sumsq4(q_aux_2);
        if usePrior
            q_aux_3 = slp*Px(W2, p, priorNorm);
            qq = qq + sumsq4(q_aux_3);
        end
        alpha=gamma/qq;

        x=x+alpha*p;

        if measurequality
            qualMeasOut(:,iter)=Measure_Quality(x_prev,x,QualMeasOpts);
        end

        % The explicit residual can drift from the recursive one CGLS keeps,
        % through the mismatch of the backprojection w.r.t. the real adjoint
        % and single precision arithmetic.
        if computeResidual
            aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
            resL2(iter)=im3Dnorm(aux,'L2');
        end
        if checkResidual
            fullres=double(resL2(iter))^2+double(sumsq4(sl*Lx(W1,x)));
            if usePrior
                fullres=fullres+double(sumsq4(slp*Px(W2,x-prior,priorNorm)));
            end
            % Rises at the level of single precision rounding are not a loss
            % of orthogonality; undoing those steps only throws work away.
            if fullres>fullres_prev*(1+1e-4)
                % we lost orthogonality: undo the bad step
                x=x-alpha*p;
                iter=iter-1;
                if restart && restarted_at~=iter
                    % restart CGLS from the current image, same weights
                    restarted_at=iter;
                    restartCGLS=true;
                    if verbose
                        disp(['Orthogonality lost, restarting CGLS at iteration ', num2str(iter)])
                    end
                    continue
                end
                % restarting did not help (or was not allowed): reweight
                if verbose
                    disp(['Orthogonality lost at iteration ', num2str(iter), ', reweighting'])
                end
                break
            end
            fullres_prev=fullres;
        end

        % If step is adequate, then continue with CGLS
        r_aux_1 = r_aux_1-alpha*q_aux_1;
        r_aux_2 = r_aux_2-alpha*q_aux_2;
        s = Atb(r_aux_1 ,geo,angles,'matched','gpuids',gpuids) + sl*Ltx(W1, r_aux_2);
        if usePrior
            r_aux_3 = r_aux_3-alpha*q_aux_3;
            s = s + slp*Ptx(W2, r_aux_3, priorNorm);
        end

        gamma1=norm(s(:),2)^2;
        beta=gamma1/gamma;
        gamma=gamma1;
        p=s+beta*p;

        % replaces double explicit loops. Breaks the inner loop when its
        % time to do a cold restart.
        if mod(iter,niter_break)==0
            break
        end
        if (iter==1 && verbose)
            expected_time=toc*niter;
            disp(name);
            disp(['Expected duration   :    ',secs2hms(expected_time)]);
            disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
            disp('');
        end
    end
    x_outer=x;
    if iter==iter_outer_start
        % a whole outer iteration without a single accepted step: more
        % reweighting from the same image will not change anything
        if verbose
            disp(['Algorithm stopped in iteration ', num2str(iter), ': no further progress.'])
        end
        break
    end


end
end

function s = sumsq4(q)
% squared L2 norm, accumulated per volume as the rest of TIGRE does
s=0;
for k=1:size(q,4)
    s=s+im3Dnorm(q(:,:,:,k),'L2')^2;
end
end

function out = Px(W, x, priorNorm)
% the weighted operator of the prior term
switch priorNorm
    case 'tv'
        out = Lx(W, x);
    case 'l1'
        out = W .* x;
    otherwise % 'l2'
        out = x;
end
end

function out = Ptx(W, y, priorNorm)
% its transpose
switch priorNorm
    case 'tv'
        out = Ltx(W, y);
    case 'l1'
        out = W .* y;
    otherwise % 'l2'
        out = y;
end
end

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

function out = Lx (W, x)
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
end

function out = Ltx (W, x)
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
end
%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids,lambda,lambda_prior,priorNorm,niter_outer,gt,restart,checkResidual,outerRestart,W1,W2]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','lambda','lambda_prior','priornorm','niter_outer','groundtruth','restart','checkresidual','outerrestart','weights','weightsprior'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:IRN_PICCS_CGLS:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:IRN_PICCS_CGLS:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
                error('TIGRE:IRN_PICCS_CGLS:InvalidInput','Invalid Init option')
            end
        case 'niter_outer'
            if default
                niter_outer=10;
            else
                niter_outer=val;
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
                    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','Invalid image for initialization');
                end
            end
            %  =========================================================================
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE:Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
        case 'lambda'
            if default
                lambda=10;
            else
                lambda=val;
            end
        case 'groundtruth'
            if default
                gt=nan;
            else
                gt=val;
            end
        case 'restart'
            if default
                restart=true;
            else
                restart=val;
            end
        case 'checkresidual'
            if default
                checkResidual=true;
            else
                checkResidual=val;
            end
        case 'lambda_prior'
            if default
                lambda_prior=10;
            else
                lambda_prior=val;
            end
        case 'priornorm'
            if default
                priorNorm='tv';
            else
                priorNorm=lower(val);
                if ~any(strcmp(priorNorm,{'tv','l2','l1'}))
                    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','priorNorm should be ''tv'', ''l2'' or ''l1''');
                end
            end
        case 'outerrestart'
            if default
                outerRestart='cold';
            else
                outerRestart=lower(val);
                if ~any(strcmp(outerRestart,{'cold','warm'}))
                    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','outerRestart should be ''cold'' or ''warm''');
                end
            end
        case 'weights'
            if default
                W1=[];
            else
                if ~isequal(size(val),geo.nVoxel')
                    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','weights should be of size geo.nVoxel''');
                end
                W1=single(val);
            end
        case 'weightsprior'
            if default
                W2=[];
            else
                if ~isequal(size(val),geo.nVoxel')
                    error('TIGRE:IRN_PICCS_CGLS:InvalidInput','weightsPrior should be of size geo.nVoxel''');
                end
                W2=single(val);
            end
        otherwise
            error('TIGRE:IRN_PICCS_CGLS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in IRN_PICCS_CGLS() or IRN_TV_CGLS()']);
    end
end


end