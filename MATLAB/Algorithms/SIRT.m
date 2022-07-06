function [res,errorL2,qualMeasOut]=SIRT(proj,geo,angles,niter,varargin)
% SIRT solves Cone Beam CT image reconstruction using Oriented Subsets
%              Simultaneous Algebraic Reconxtruction Techique algorithm
%
%   SIRT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%   SIRT(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%
%   'lambda':      Sets the value of the hyperparameter. Default is 1
%
%   'lambda_red':   Reduction of lambda.Every iteration
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
%   'InitImg'      an image for the 'image' initialization. Avoid.
%
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%   'QualMeas'     Asks the algorithm for a set of quality measurement
%                  parameters. Input should contain a cell array of desired
%                  quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                  These will be computed in each iteration.
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

%% Deal with input parameters

[lambda,res,lambdared,verbose,QualMeasOpts,nonneg,gpuids]=parse_inputs(proj,geo,angles,varargin);
measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter);

if nargout>1
    computeL2=true;
else
    computeL2=false;
end
errorL2=[];

%% initialize stuff

%% Create weigthing matrices

% Projection weigth, W

geoaux=geo;
geoaux.sVoxel([1 2])=geo.sVoxel([1 2])*1.1; % a Bit bigger, to avoid numerical division by zero (small number)
geoaux.sVoxel(3)=max(geo.sDetector(2),geo.sVoxel(3)); % make sure lines are not cropped. One is for when image is bigger than detector and viceversa
geoaux.nVoxel=[2,2,2]'; % accurate enough?
geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'Siddon','gpuids',gpuids);  %
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;
clear geoaux;

% Back-Projection weigth, V
 V=computeV(geo,angles,{angles},{1:length(angles)},'gpuids',gpuids);

%% hyperparameter stuff
nesterov=false;
if ischar(lambda)&&strcmp(lambda,'nesterov')
    nesterov=true;
    lambda=(1+sqrt(1+4))/2;
    gamma=0;
    ynesterov=zeros(size(res),'single');
    ynesterov_prev=ynesterov;
end
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
    % --------- Memory expensive-----------
    %         proj_err=proj-Ax(res,geo,angles);                 %                                 (b-Ax)
    %         weighted_err=W.*proj_err;                         %                          W^-1 * (b-Ax)
    %         backprj=Atb(weighted_err,geo,angles);             %                     At * W^-1 * (b-Ax)
    %         weigth_backprj=bsxfun(@times,1./V,backprj);       %                 V * At * W^-1 * (b-Ax)
    %         res=res+lambda*weigth_backprj;                    % x= x + lambda * V * At * W^-1 * (b-Ax)
    % ------------------------------------
    % --------- Memory cheap(er)-----------
    
    if nesterov
        % The nesterov update is quite similar to the normal update, it
        % just uses this update, plus part of the last one.
        ynesterov=res + bsxfun(@times,1./V,Atb(W.*(proj-Ax(res,geo,angles,'gpuids',gpuids)),geo,angles,'gpuids',gpuids));
        res=(1-gamma)*ynesterov+gamma*ynesterov_prev;
    else
        res=res+lambda*bsxfun(@times,1./V,Atb(W.*(proj-Ax(res,geo,angles,'gpuids',gpuids)),geo,angles,'gpuids',gpuids)); % x= x + lambda * V * At * W^-1 * (b-Ax)
    end
    % ------------------------------------
    
    if nonneg
        res=max(res,0);
    end
    if nesterov
        gamma=(1-lambda);
        lambda=(1+sqrt(1+4*lambda^2))/2;
        gamma=gamma/lambda;
    else
        lambda=lambda*lambdared;
    end
    % If quality is being measured
    if measurequality
        % HERE GOES
        qualMeasOut(:,ii)=Measure_Quality(res_prev,res,QualMeasOpts);
    end
    
    if computeL2 || nesterov
        errornow=im3Dnorm(proj-Ax(res,geo,angles,'gpuids',gpuids),'L2','gpuids',gpuids);                       % Compute error norm2 of b-Ax
        % If the error is not minimized.
        if  ii~=1 && errornow>errorL2(end)
            if verbose
                disp(['Convergence criteria met, exiting on iteration number:', num2str(ii)]);
            end
            return;
        end
        errorL2=[errorL2 errornow];
    end
    
    
    if (ii==1 && verbose==1);
        expected_time=toc*niter;
        disp('SIRT');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
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
    initres=SIRT(proj,geo,alpha,niter,'Init','image','InitImg',initres,'Verbose',0,'gpuids',gpuids);
    
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


function [lambda,res,lambdared,verbose,QualMeasOpts,nonneg,gpuids]=parse_inputs(proj,geo,alpha,argin)
opts=     {'lambda','init','initimg','verbose','lambda_red','qualmeas','nonneg','gpuids'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:SIRT:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:SIRT:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:SIRT:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
        end
        val=argin{jj};
    end
    
    switch opt
        % % % % % % % Verbose
        case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE: Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
            % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=1;
            elseif ischar(val)&&strcmpi(val,'nesterov');
                lambda='nesterov'; %just for lowercase/upercase
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:SIRT:InvalidInput','Invalid lambda')
            else
                lambda=val;
            end
        case 'lambda_red'
            if default
                lambdared=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:SIRT:InvalidInput','Invalid lambda')
                end
                lambdared=val;
            end
        case 'init'
            res=[];
            if default || strcmp(val,'none')
                res=zeros(geo.nVoxel','single');
                continue;
            end
            if strcmp(val,'FDK')
                res=FDK(proj,geo,alpha);
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
                error('TIGRE:SIRT:InvalidInput','Invalid Init option')
            end
            % % % % % % % ERROR
        case 'initimg'
            if default
                continue;
            end
            if exist('initwithimage','var');
                if isequal(size(val),geo.nVoxel');
                    res=single(val);
                else
                    error('TIGRE:SIRT:InvalidInput','Invalid image for initialization');
                end
            end
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:SIRT:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'nonneg'
            if default
                nonneg=true;
            else
                nonneg=val;
            end
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
        otherwise
            error('TIGRE:SIRT:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in SIRT()']);
    end
end

end