function [x,resL2,qualMeasOut]= LSMR(proj,geo,angles,niter,varargin)

% LSMR solves the CBCT problem using LSMR.
%
%  LSMR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%  LSMR(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%  'lambda'  Value of parameter lambda, default 0.
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
%  'restart'  true or false. By default the algorithm will restart when
%             loss of ortogonality is found.
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

[verbose,x,QualMeasOpts,gpuids,lambda,gt,restart]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    x0=gt;
    clear gt
end
if nargout<3 && measurequality
    warning("Image metrics requested but none catched as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),niter);
resL2 = zeros(1,niter);

% David Chin-Lung Fong and Michael Saunders //doi.org/10.1137/10079687X


% Enumeration as given in the paper for 'Algorithm LSMR'
iter=0;
remember=0;
while iter<niter
    % (1) Initialize
    u=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
    normr = norm(u(:),2);
    beta = normr;
    u = u/beta;
    
    
    v=Atb(u,geo,angles,'matched','gpuids',gpuids);
    alpha = norm(v(:),2); % msl: do we want to check if it is 0?
    v = v/alpha;
    
    alphabar = alpha;
    zetabar = alpha * beta;
    rho = 1;
    rhobar = 1;
    cbar = 1;
    sbar = 0;
    h = v;
    hbar = 0;
    
    % Compute the residual norm ||r_k||
    betadd = beta;
    betad = 0;
    rhod = 1;
    tautilda = 0;
    thetatilda = 0;
    zeta = 0;
    d = 0;
    
    % (2) Start iterations
    for ii=iter:niter
        iter=iter+1;
        if measurequality && ~strcmp(QualMeasOpts,'error_norm')
            x0 = x; % only store if necessary
        end
        if (iter==1 && verbose);tic;end
        
        % (3) Continue the bidiagonalization
        u = Ax(v,geo,angles,'Siddon','gpuids',gpuids) - alpha*u;
        beta = norm(u(:),2);
        u = u / beta;
        
        v = Atb(u,geo,angles,'matched','gpuids',gpuids) - beta*v;
        
        alpha = norm(v(:),2);
        v = v / alpha;
        
        % (4) Construct and apply rotation \hat{P}_k
        alphahat = sqrt(alphabar^2 + lambda^2);
        chat = alphabar/alphahat;
        shat = lambda/alphahat;
        
        % (5) Construct and apply rotation P_k
        rhopre = rho;
        rho = sqrt(alphahat^2 + beta^2);
        c = alphahat / rho;
        s = beta / rho;
        theta = s * alpha;
        alphabar = c * alpha;
        
        % (6) Construct and apply rotation \bar{P}_k
        thetabar = sbar * rho;
        rhobarpre = rhobar;
        rhobar = sqrt((cbar*rho)^2 + theta^2);
        cbar = cbar * rho / rhobar;
        sbar = theta / rhobar;
        zetapre = zeta;
        zeta = cbar * zetabar;
        zetabar = -sbar * zetabar;
        
        % (7) Update \bar{h}, x, h
        hbar = h - (thetabar*rho/(rhopre*rhobarpre))*hbar;
        x = x + (zeta / (rho*rhobar)) * hbar;
        h = v - (theta / rho) * h;
        
        % (8) Apply rotation \hat{P}_k, P_k
        betaacute = chat* betadd;
        betacheck = - shat* betadd;
        
        % Computing ||r_k||
        
        betahat = c * betaacute;
        betadd = -s * betaacute;
        
        % Update estimated quantities of interest.
        %  (9) If k >= 2, construct and apply \tilde{P} to compute ||r_k||
        rhotilda = sqrt(rhod^2 + thetabar^2);
        ctilda = rhod / rhotilda;
        stilda = thetabar / rhotilda;
        thetatildapre = thetatilda;
        thetatilda = stilda * rhobar;
        rhod = ctilda * rhobar;
        % betatilda = ctilda * betad + stilda * betahat; % msl: in the orinal paper, but not used
        betad = -stilda * betad + ctilda * betahat;
        
        % (10) Update \tilde{t}_k by forward substitution
        tautilda = (zetapre - thetatildapre*tautilda) / rhotilda;
        taud = (zeta - thetatilda*tautilda) / rhod;
        
        % (11) Compute ||r_k||
        d = d + betacheck^2;
        gamma_var = d + (betad - taud)^2 + betadd^2;
        aux = sqrt(gamma_var); % this is the residual teh algorithm follows, but we lose ortogonality, so we compute it explicitly
        
        % ||A^T r_k || is just |zetabar|
        
        
        
        % (6) Test for convergence.
        % msl: I still need to implement this.
        % msl: There are suggestions on the original paper. Let's talk about it!
        
        if measurequality
            qualMeasOut(:,iter)=Measure_Quality(x0,x,QualMeasOpts);
        end
        % The following should never happen, but the reallity is that if we use
        % the residual from the algorithm, it starts diverging from this explicit residual value.
        % This is an interesting fact that I believe may be caused either by
        % the mismatch of the backprojection w.r.t the real adjoint, or
        % numerical issues related to doing several order of magnitude
        % difference operations on single precission numbers.
        aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
        resL2(iter)=im3Dnorm(aux,'L2');
        if iter>1 && resL2(iter)>resL2(iter-1)
            % we lost orthogonality, lets restart the algorithm unless the
            % user asked us not to.
            
            % undo bad step.
            x=x-(zeta / (rho*rhobar)) * hbar;
            % if the restart didn't work.
            if remember==iter || ~restart
                disp(['Algorithm stoped in iteration ', num2str(iter),' due to loss of ortogonality.'])
                return;
            end
            remember=iter;
            iter=iter-1;
            if verbose
                disp(['Orthogonality lost, restarting at iteration ', num2str(iter) ])
            end
            break
        end
        
        if (iter==1 && verbose)
            expected_time=toc*niter;
            disp('LSMR');
            disp(['Expected duration   :    ',secs2hms(expected_time)]);
            disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
            disp('');
        end
    end
end
end

%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids, lambda,gt,restart]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','lambda','groundtruth','restart'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:LSMR:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:LSMR:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:LSMR:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
                error('TIGRE:LSMR:InvalidInput','Invalid Init option')
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
                    error('TIGRE:LSMR:InvalidInput','Invalid image for initialization');
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
                    error('TIGRE:LSMR:InvalidInput','Invalid quality measurement parameters');
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
                lambda = 0;
            else
                lambda = val;
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
        otherwise
            error('TIGRE:LSMR:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in LSMR()']);
    end
end


end
