function [x,residual,qualMeasOut]= LSMR(proj,geo,angles,niter,varargin)

% LSMR solves the CBCT problem using LSMR.
% 
%  LSMR(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
% 
%  LSMR(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
% 
%  'lambda'  Value of parameter lambda, default 0. 
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
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Malena Sabate Landman, Ander Biguri 
%--------------------------------------------------------------------------
%%

[verbose,x,QualMeasOpts,gpuids,lambda]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter);


% David Chin-Lung Fong and Michael Saunders //doi.org/10.1137/10079687X

% Enumeration as given in the paper for 'Algorithm LSMR'
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

% msl: is this error the residual norm ? 
residual = zeros(1,niter); 

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
    residual(ii) = sqrt(gamma_var);
    
    % ||A^T r_k || is just |zetabar|


  
    % (6) Test for convergence. 
    % msl: I still need to implement this. 
    % msl: There are suggestions on the original paper. Let's talk about it!
    
    if measurequality % msl: what is this??
        qualMeasOut(:,ii)=Measure_Quality(x0,x,QualMeasOpts);
    end
     
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
function [verbose,x,QualMeasOpts,gpuids, lambda]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','lambda'};
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
        otherwise 
            error('TIGRE:LSMR:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in LSMR()']);
    end
end


end
