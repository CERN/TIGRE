function [res,resL2,qualMeasOut]=OSEM(proj,geo,angles,niter,varargin)
% OS_SART solves Cone Beam CT image reconstruction using Oriented Subsets
%              Simultaneous Algebraic Reconstruction Technique algorithm
%
%   OS_SART(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%   OS_SART(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%   'BlockSize':   Sets the projection block size used simultaneously. If
%                   if BlockSize = length(alpha) then OSEM becomes MLEM. Default is 20.
%
%   'lambda':      Sets the value of the hyperparameter. Default is 1
%
%   'lambda_red':  Reduction of lambda. Every iteration
%                  lambda=lambdared*lambda. Default is 0.95
%
%   'Init':        Describes different initialization techniques.
%                  'none'     : Initializes the image to zeros (default)
%                  'FDK'      : Initializes image to FDK reconstruction
%                  'multigrid': Initializes image by solving the problem in
%                               small scale and increasing it when relative
%                               convergence is reached.
%                  'image'    : Initialization using a user specified
%                               image. Not recommended unless you really
%                               know what you are doing.
%   'InitImg'      an image for the 'image' initialization. Avoid.
%
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%   'QualMeas'     Asks the algorithm for a set of quality measurement
%                  parameters. Input should contain a cell array of desired
%                  quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                  These will be computed in each iteration.
% 'OrderStrategy'  Chooses the subset ordering strategy. Options are
%                  'ordered' : uses them in the input order, but divided
%                  'random'  : orders them randomly
%                  'angularDistance': chooses the next subset with the
%                                     biggest angular distance with the ones used.
% 'redundancy_weighting': true or false. Default is true. Applies data
%                         redundancy weighting to projections in the update step
%                         (relevant for offset detector geometry)
%  'groundTruth'  an image as ground truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
% OUTPUTS:
%
%    [img]                       will output the reconstructed image
%    [img,errorL2]               will output the L2 norm of the residual
%                                (the function being minimized)
%    [img,errorL2,qualMeas]      will output the quality measurements asked
%                                by the input 'QualMeas'
%
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

[blocksize,lambda,res,lambdared,verbose,QualMeasOpts,OrderStrategy,nonneg,gpuids,redundancy_weights,gt]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    res_prev=gt;
    clear gt
end
if nargout<3 && measurequality
    warning("Image metrics requested but none catched as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),niter);

resL2=zeros(1,niter);


if nargout>1
    computeL2=true;
else
    computeL2=false;
end
% does detector rotation exists?
if ~isfield(geo,'rotDetector')
    geo.rotDetector=[0;0;0];
end
%% weight matrices
% first order the projection angles
[alphablocks,orig_index]=order_subsets(angles,blocksize,OrderStrategy);


% Projection weight, W
W=computeW(geo,angles,gpuids);
% Back-Projection weight, V
V=computeV(geo,angles,alphablocks,orig_index,'gpuids',gpuids);

if redundancy_weights
    % Data redundancy weighting, W_r implemented using Wang weighting
    % reference: https://iopscience.iop.org/article/10.1088/1361-6560/ac16bc
    
    num_frames = size(proj,3);
    W_r = redundancy_weighting(geo);
    W_r = repmat(W_r,[1,1,num_frames]);
    % disp('Size of redundancy weighting matrix');
    % disp(size(W_r));
    W = W.*W_r; % include redundancy weighting in W
end

clear A x y dx dz


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
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;
rotDetector=geo.rotDetector;
DSD=geo.DSD;
DSO=geo.DSO;


% TODO : Add options for Stopping criteria
for ii=1:niter
    
    % If verbose, time the algorithm
    if (ii==1 && verbose==1);tic;end
    % If quality is going to be measured, then we need to save previous image
    % THIS TAKES MEMORY!
    if measurequality && ~strcmp(QualMeasOpts,'error_norm')
        res_prev = res; % only store if necessary
    end
    
    
    for jj=1:length(alphablocks)
        % Get offsets
        if size(offOrigin,2)==size(angles,2)
            geo.offOrigin=offOrigin(:,orig_index{jj});
        end
        if size(offDetector,2)==size(angles,2)
            geo.offDetector=offDetector(:,orig_index{jj});
        end
        if size(rotDetector,2)==size(angles,2)
            geo.rotDetector=rotDetector(:,orig_index{jj});
        end
        if size(DSO,2)==size(angles,2)
            geo.DSO=DSO(:,orig_index{jj});
        end
        if size(DSD,2)==size(angles,2)
            geo.DSD=DSD(:,orig_index{jj});
        end
        den = Ax(res,geo,alphablocks{:,jj},'gpuids',gpuids);
        den(den<=0.)=inf;
        
        imgupdate = Atb(proj(:,:,orig_index{jj})./den, geo,alphablocks{:,jj},'matched','gpuids',gpuids)./sum(V(:,:,jj),3);
        
        res=res.*imgupdate;
        % Non-negativity constraint
        if nonneg
              res = max(res,0.);
        end
        
    end
    
    % reduce hyperparameter
    if nesterov
        gamma=(1-lambda);
        lambda=(1+sqrt(1+4*lambda^2))/2;
        gamma=gamma/lambda;
    else
        lambda=lambda*lambdared;
    end
    
    if measurequality
        qualMeasOut(:,ii)=Measure_Quality(res_prev,res,QualMeasOpts);
    end
    
    if computeL2 || nesterov
        % Compute error norm2 of b-Ax
        geo.offOrigin=offOrigin;
        geo.offDetector=offDetector;
        geo.DSD=DSD;
        geo.rotDetector=rotDetector;
        resL2(ii)=im3Dnorm(proj-Ax(res,geo,angles,'Siddon','gpuids',gpuids),'L2');
        % If the error is not minimized
        if ii~=1 && resL2(ii)>resL2(ii-1)
            if verbose
                disp(['Convergence criteria met, exiting on iteration number:', num2str(ii)]);
            end
            return
        end
        % Store Error
        
    end
    % If timing was asked for
    if ii==1 && verbose==1
        expected_time=toc*(niter-1);
        expected_duration=toc*(niter);
        disp('OSEM');
        disp(['Expected duration   :    ',secs2hms(expected_duration)]);
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
    return
end
niter=100;
nblock=20;
initres=zeros(geo.nVoxel');
while ~isequal(geo.nVoxel,finalsize)
    
    
    % solve subsampled grid
    initres=OS_SART(proj,geo,alpha,niter,'BlockSize',nblock,'Init','image','InitImg',initres,'Verbose',0);
    
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

%% Parse inputs
function [block_size,lambda,res,lambdared,verbose,QualMeasOpts,OrderStrategy,nonneg,gpuids,redundancy_weights,gt]=parse_inputs(proj,geo,alpha,argin)
opts={'blocksize','lambda','init','initimg','verbose','lambda_red','qualmeas','orderstrategy','nonneg','gpuids','redundancy_weighting','groundtruth'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:OS_SART:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:OS_SART:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option is not default, then extract value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end
        if isempty(ind)
            error('TIGRE:OS_SART:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
            else
                lambda=val;
            end
        case 'lambda_red'
            if default
                lambdared=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
                end
                lambdared=val;
            end
        case 'blocksize'
            if default
                block_size=20;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:OS_SART:InvalidInput','Invalid BlockSize')
                end
                block_size=val;
            end
            
        case 'init'
            res=[];
            if default || strcmp(val,'none')
                res=ones(geo.nVoxel','single');
                continue
            end
            if strcmp(val,'FDK')
                res=FDK(proj,geo,alpha);
                continue
            end
            if strcmp(val,'multigrid')
                res=init_multigrid(proj,geo,alpha);
                continue
            end
            if strcmp(val,'image')
                initwithimage=1;
                continue
            end
            if isempty(res)
                error('TIGRE:OS_SART:InvalidInput','Invalid Init option')
            end
            % % % % % % % ERROR
        case 'initimg'
            if default
                continue
            end
            if exist('initwithimage','var')
                if isequal(size(val),geo.nVoxel')
                    res=single(val);
                else
                    error('TIGRE:OS_SART:InvalidInput','Invalid image for initialization');
                end
            end
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:OS_SART:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'orderstrategy'
            if default
                OrderStrategy='random';
            else
                OrderStrategy=val;
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
        case 'redundancy_weighting'
            if default
                redundancy_weights = true;
            else
                redundancy_weights = val;
            end
        case 'groundtruth'
            if default
                gt=nan;
            else
                gt=val;
            end
        otherwise
            error('TIGRE:OS_SART:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in OS_SART_CBCT()']);
    end
end

end