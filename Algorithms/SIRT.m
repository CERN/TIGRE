function [res,errorL2,QualMeasOpts]=SIRT(proj,geo,alpha,niter,varargin)
% SIRT_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets
%              Simultaneous Algebraic Reconxtruction Techique algorithm
%
%   SIRT_CBCT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%   SIRT_CBCT(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%
%   'lambda':      Sets the value of the hyperparameter. Default is 1
%
%   'lambdared':   Reduction of lambda.Every iteration
%                  lambda=lambdared*lambda. Default is 0.95
%
%   'Init':        Describes diferent initialization techniques.
%                  'none'     : Initializes the image to zeros (default)
%                  'FDK'      : intializes image to FDK reconstrucition
%                  'multigrid': Initializes image by solving the problem in
%                               small scale and increasing it when relative
%                               convergence is reached.
%                  'image'    : Initialization using a user specified
%                               image. Not recomended unless you really
%                               know what you are doing.
%   'InitImg'      an image for the 'image' initialization. Aviod.
%
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.

%% Deal with input parameters

[lambda,res,lamdbared,verbose,QualMeasOpts]=parse_inputs(proj,geo,alpha,varargin);
measurequality=~isempty(QualMeasOpts);
errorL2=[];

%% initialize stuff

%% Create weigthing matrices

% Projection weigth, W
% Projection weigth, W
W=Ax(ones(geo.nVoxel'),geo,alpha,'Krylov');  %
W(W<min(geo.dVoxel))=Inf;
W=1./W;
% Back-Projection weigth, V
[x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
    -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
A = permute(alpha, [1 3 2]);
V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
V=sum(V,3);
clear A x y dx dz;

%% Iterate

errorL2=[];

% TODO : Add options for Stopping criteria
for ii=1:niter
    if (ii==1 && verbose==1);tic;end
    % If quality is going to be measured, then we need to save previous image
    % THIS TAKES MEMORY!
    if measurequality
        res_prev=res;
    end
       
    proj_err=proj-Ax(res,geo,alpha);                  %                                 (b-Ax)
    weighted_err=W.*proj_err;                         %                          W^-1 * (b-Ax)
    backprj=Atb(weighted_err,geo,alpha);              %                     At * W^-1 * (b-Ax)
    weigth_backprj=bsxfun(@times,1./V,backprj);       %                 V * At * W^-1 * (b-Ax)
    res=res+lambda*weigth_backprj;                    % x= x + lambda * V * At * W^-1 * (b-Ax)
    res(res<0)=0;
    
   % If quality is being measured
    if measurequality
       % HERE GOES  
       qualMeas=Measure_Quality(res_prev,res,QualMeasOpts);
    end   
    
    errornow=norm(proj_err(:));                       % Compute error norm2 of b-Ax
    % If the error is not minimized.
    if ii>1 && errornow>errorL2(end)
        return;
    end
    lambda=lambda*lamdbared;
    
    errorL2=[errorL2 errornow];
    if (ii==1 && verbose==1);
        expected_time=toc*niter;
        disp('SIRT');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
end





end


function initres=init_multigrid(proj,geo,alpha)

finalsize=geo.nVoxel;
% start with 64
geo.nVoxel=[64;64;64];
geo.dVoxel=geo.sVoxel./geo.nVoxel;
if any(finalsize<geo.nVoxel)
    initres=zeros(finalsize');
    return;
end
niter=100;
initres=zeros(geo.nVoxel');
while ~isequal(geo.nVoxel,finalsize)
    
    
    % solve subsampled grid
    initres=SIRT_CBCT(proj,geo,alpha,niter,'Init','image','InitImg',initres,'Verbose',0);
    
    % Get new dims.
    geo.nVoxel=geo.nVoxel*2;
    geo.nVoxel(geo.nVoxel>finalsize)=finalsize(geo.nVoxel>finalsize);
    geo.dVoxel=geo.sVoxel./geo.nVoxel;
    % Upsample!
    % (hopefully computer has enough memory............)
    [y, x, z]=ndgrid(linspace(1,size(initres,1),geo.nVoxel(1)),...
        linspace(1,size(initres,2),geo.nVoxel(2)),...
        linspace(1,size(initres,3),geo.nVoxel(3)));
    initres=interp3(initres,x,y,z);
    clear x y z
end


end


function [lambda,res,lamdbared,verbose,QualMeasOpts]=parse_inputs(proj,geo,alpha,argin)
opts=     {'lambda','Init','InitImg','Verbose','lambdaRed','QualMeas'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('CBCT:SIRT:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,argin{ii}));
    if ~isempty(ind)
        defaults(ind)=0;
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,argin{jj}));
            jj=jj+1;
        end
        val=argin{jj};
    end
    
    switch opt
        % % % % % % % Verbose
        case 'Verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=0.95;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:SIRT:InvalidInput','Invalid lambda')
                end
                lambda=val;
            end
        case 'lambdaRed'
            if default
                lamdbared=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:SIRT:InvalidInput','Invalid lambda')
                end
                lamdbared=val;
            end
        case 'Init'
            res=[];
            if default || strcmp(val,'none')
                res=zeros(geo.nVoxel');
                continue;
            end
            if strcmp(val,'FDK')
                res=FDK_CBCT(proj,geo,alpha);
                continue;
            end
            if strcmp(val,'multigrid')
                res=init_multigrid(proj,geo,alpha);
                continue;
            end
            if strcmp(val,'image')
                initwithimage=1;
                continue;
            end
            if isempty(res)
                error('CBCT:SIRT:InvalidInput','Invalid Init option')
            end
            % % % % % % % ERROR
        case 'InitImg'
            if default
                continue;
            end
            if exist('initwithimage','var');
                if isequal(size(val),geo.nVoxel');
                    res=val;
                else
                    error('CBCT:SIRT:InvalidInput','Invalid image for initialization');
                end
            end
        case 'QualMeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('CBCT:SIRT:InvalidInput','Invalid quality measurement parameters');
                end
            end
        otherwise
            error('CBCT:SIRT:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in SIRT()']);
    end
end

end