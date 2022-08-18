function [x,errorL2,qualMeasOut]= LSQR(proj,geo,angles,niter,varargin)

% LSQR_CBCT solves the CBCT problem using the conjugate gradient least
% squares. This is mathematically equivalent to CGLS_CBCT.
% 
%  LSQR_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
% 
%  LSQR_CBCT(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
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
% 'redundancy_weighting': true or false. Default is true. Applies data
%                         redundancy weighting to projections in the update step
%                         (relevant for offset detector geometry)
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
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------

%%

% msl: I assume this x is either an array of 0 of the appropriate
% dimensions or an initial guess

[verbose,x,QualMeasOpts,gpuids,redundancy_weights]=parse_inputs(proj,geo,angles,varargin);

% msl: no idea of what this is. Should I check?
measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter);
if redundancy_weights
    % Data redundancy weighting, W_r implemented using Wang weighting
    % reference: https://iopscience.iop.org/article/10.1088/1361-6560/ac16bc
    
    num_frames = size(proj,3);
    W_r = redundancy_weighting(geo);
    W_r = repmat(W_r,[1,1,num_frames]);
    % disp('Size of redundancy weighting matrix');
    % disp(size(W_r));
end

% Paige and Saunders //doi.org/10.1145/355984.355989

% Enumeration as given in the paper for 'Algorithm LSQR'
% (1) Initialize 
u=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
normr = norm(u(:),2);
u = u/normr;
% msl: why is r not notmalized in CGLS? Check.

beta = normr;
phibar = beta;

% msl: I have not checked 'redundancy_weights', do they have to do with 
% mismatched transpose? 
if redundancy_weights 
    v=Atb(W_r .* u,geo,angles,'matched','gpuids',gpuids);
else
    v=Atb(u,geo,angles,'matched','gpuids',gpuids);
end

alpha = norm(v(:),2); % msl: do we want to check if it is 0? 
v = v/alpha;
rhobar = alpha;
w = v;

normAtr = beta*alpha; % Do we want this? ||A^T r_k||
% msl: do we want to check for convergence? In well posed problems it would
% make sense, not sure now.

errorL2=zeros(1,niter); % msl: is this error the residual norm ? 

% (2) Start iterations 
for ii=1:niter
    x0 = x;
    if (ii==1 && verbose);tic;end
    
    % (3)(a)
    u = Ax(v,geo,angles,'Siddon','gpuids',gpuids) - alpha*u;
    beta = norm(u(:),2);
    u = u / beta;
    
    % (3)(b)
    if redundancy_weights % msl: Is this right here?
        v = Atb(W_r .* u,geo,angles,'matched','gpuids',gpuids) - beta*v;
    else
        v = Atb(u,geo,angles,'matched','gpuids',gpuids) - beta*v;
    end
    alpha = norm(v(:),2);
    v = v / alpha;    

    % (4)(a-g)
    rho = sqrt(rhobar^2 + beta^2);
    c = rhobar / rho;
    s = beta / rho;
    theta = s * alpha;
    rhobar = - c * alpha; % msl: I am overwriting a matlab function, maybe change name
    phi = c * phibar;
    phibar = s * phibar;
    
    % (5) Update x, w
    x = x + (phi / rho) * w;
    w = v - (theta / rho) * w;

    % Update estimated quantities of interest.
    % msl: We can also compute cheaply estimates of ||x||, ||A||, cond(A)
    normr = normr*abs(s);   % ||r_k|| = ||b - A x_k||
                            % Only exact if we do not have orth. loss 
    normAtr = phibar * alpha * abs(c); % msl: Do we want this? ||A^T r_k||
    errorL2(ii)=normr;
    
    % (6) Test for convergence. 
    % msl: I still need to implement this. 
    % msl: There are suggestions on the original paper. Let's talk about it!
    
    if measurequality % msl: what is this??
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
function [verbose,x,QualMeasOpts,gpuids,redundancy_weights]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','redundancy_weighting'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:CGLS:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
       error('TIGRE:CGLS:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]); 
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
            error('TIGRE:CGLS:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
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
               error('TIGRE:CGLS:InvalidInput','Invalid Init option') 
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
                    error('TIGRE:CGLS:InvalidInput','Invalid image for initialization');
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
                    error('TIGRE:CGLS:InvalidInput','Invalid quality measurement parameters');
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
        case 'redundancy_weighting'
            if default
                redundancy_weights = true;
            else
                redundancy_weights = val;
            end
        otherwise 
            error('TIGRE:CGLS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end