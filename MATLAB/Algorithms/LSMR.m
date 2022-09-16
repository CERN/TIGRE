function [x,errorL2,qualMeasOut]= LSMR(proj,geo,angles,niter,varargin)

% LSMR_CBCT solves the CBCT problem using LSMR.
% 
%  LSMR_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
% 
%  LSMR_CBCT(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
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
% Coded by:           Malena Sabate Landman, Ander Biguri 
%--------------------------------------------------------------------------
%%

[verbose,x,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter);


% David Chin-Lung Fong and Michael Saunders //doi.org/10.1137/10079687X

% Enumeration as given in the paper for 'Algorithm LSMR'
% (1) Initialize 
u=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
normr = norm(u(:),2);
u = u/normr;

beta = normr;

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


% msl: is this error the residual norm ? 
errorL2=zeros(1,niter); 

% (2) Start iterations 
for ii=1:niter
    x0 = x;
    if (ii==1 && verbose);tic;end
    
    % (3) Continue the bidiagonalization
    u = Ax(v,geo,angles,'Siddon','gpuids',gpuids) - alpha*u;
    beta = norm(u(:),2);
    u = u / beta;
    
    v = Atb(u,geo,angles,'matched','gpuids',gpuids) - beta*v;

    alpha = norm(v(:),2);
    v = v / alpha;    

    % (4) Construct and apply rotation Pk
    rhopre = rho;
    rho = sqrt(alphabar^2 + beta^2);
    c = alphabar / rho;
    s = beta / rho;
    theta = s * alpha;
    alphabar = c * alpha;
  
    % Computing ||r_k||, eqn (3.4) in paper
    betahat = c * betadd; 
    betadd = -s * betadd;

    % (5) Construct and apply rotation Pkbar
    thetabar = sbar * rho;
    rhobarpre = rhobar;
    rhobar = sqrt((cbar*rho)^2 + theta^2);
    cbar = cbar * rho / rhobar;
    sbar = theta / rhobar;
    zetapre = zeta;
    zeta = cbar * zetabar; 
    zetabar = -sbar * zetabar;

    % (6) Update h, hbar, x
    hbar = h - (thetabar*rho/(rhopre*rhobarpre))*hbar;
    x = x + (zeta / (rho*rhobar)) * hbar;
    h = v - (theta / rho) * h;

    % Update estimated quantities of interest.
    % If k >= 2, construct and apply Ptilda to compute ||r_k||
    if ii >= 1
        % eqn (3.5) in paper
        rhotilda = sqrt(rhod^2 + thetabar^2);
        ctilda = rhod / rhotilda;
        stilda = thetabar / rhotilda;
        thetatildapre = thetatilda;
        thetatilda = stilda * rhobar;
        rhod = ctilda * rhobar;
        betatilda = ctilda * betad + stilda * betahat; % msl: in the orinal paper, but not used
        betad = -stilda * betad + ctilda * betahat;
        % Update ttilda by forward substitution
        tautilda = (zetapre - thetatildapre*tautilda) / rhotilda;
        taud = (zeta - thetatilda*tautilda) / rhod;
        % Form ||r_k||
        gamma_var = (betad - taud)^2 + betadd^2;
        errorL2(ii) = sqrt(gamma_var);
        % ||A^T r_k || is just |zetabar|
    else     
        % Not sure what to say about errorL2(1) ... 
        errorL2(ii)=NaN;
    end

    aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids); %expensive, is there any way to check this better?
    errorL2(ii)=im3Dnorm(aux,'L2');
    
    
    % (6) Test for convergence. 
    % msl: I still need to implement this. 
    % msl: There are suggestions on the original paper. Let's talk about it!
    
    if measurequality % msl: what is this??
        qualMeasOut(:,ii)=Measure_Quality(x0,x,QualMeasOpts);
    end

%   Does not make sense here, we are minimizing ||A^T rk|| not ||r_k||  
%     if ii>1 && errorL2(ii)>errorL2(ii-1)  % msl: not checked
%         % OUT!
%        x=x-alpha*v;
%        if verbose
%           disp(['CGLS stoped in iteration N', num2str(ii),' due to divergence.']) 
%        end
%        return; 
%     end
     
    if (ii==1 && verbose)
        expected_time=toc*niter;   
        disp('LSMR');
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
