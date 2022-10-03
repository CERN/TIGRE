function [x,errorL2,qualMeasOut]= hybrid_LSQR(proj,geo,angles,niter,varargin)

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

[verbose,x0,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,varargin);

% msl: no idea of what this is. Should I check?
measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter);

% Paige and Saunders //doi.org/10.1145/355984.355989

% This should be one to avoid loosing orthogonality, but can be switched
% off for testing.
reorth = 1;

% Initialise matrices
U = zeros(prod(size(proj)), niter+1);
V = zeros(prod(geo.nVoxel), niter); % Malena: Check if prod(geo.nVoxel) is correct
B = zeros(niter+1,niter); % Projected matrix
proj_rhs = zeros(niter+1,1); % Projected right hand side

% Enumeration as given in the paper for 'Algorithm LSQR'
% (1) Initialize 
u = proj-Ax(x0,geo,angles,'Siddon','gpuids',gpuids);
normr = norm(u(:),2);
u = u/normr;
U(:,1) = u(:);

beta = normr;
proj_rhs(1) = normr;
errorL2 = zeros(1,niter);

% (2) Start iterations 
for ii=1:niter

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

    % Malena: add regularization parameter choices
    lambda = 10;
    
    Dk = Sk.^2 + lambda^2;
    rhskhat = Sk .* rhskhat(1:ii);
    yhat = rhskhat(1:ii)./Dk;
    y = Vk * yhat;

    errorL2(ii)=norm(rhsk - Bk*y); % residual norm
    % norm_proj = norm (proj(:),2);
    % errorL2(ii) = norm(rhsk - Bk*y)/norm_proj; % relative residual norm

    d = V(:,1:ii)*y;
    x = x0 + reshape(d,size(x0));


    % Test for convergence. 
    % msl: I still need to implement this. 
    % msl: There are suggestions on the original paper. Let's talk about it!
    
    if measurequality 
        qualMeasOut(:,ii)=Measure_Quality(x0,x,QualMeasOpts);
    end

    if ii>1 && errorL2(ii)>errorL2(ii-1)  % msl: not checked
        % OUT!
       x=x-alpha*v;
       if verbose
          disp(['CGLS stoped in iteration N', num2str(ii),' due to divergence.']) 
       end
       return; 
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


%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids'};
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