function [res,errorL2,qualMeasOut]=OS_SART(proj,geo,angles,niter,varargin)

% OS_SART_CBCT solves Cone Beam CT image reconstruction using Oriented Subsets
%              Simultaneous Algebraic Reconxtruction Techique algorithm
%
%   OS_SART_CBCT(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%   OS_SART_CBCT(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%   'BlockSize':   Sets the projection block size used simultaneously. If
%                  BlockSize = 1 OS-SART becomes SART and if  BlockSize = length(alpha)
%                  then OS-SART becomes SIRT. Default is 20.
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
%   'QualMeas'     Asks the algorithm for a set of quality measurement
%                  parameters. Input should contain a cell array of desired
%                  quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                  These will be computed in each iteration.
% 'OrderStrategy'  Chooses the subset ordering strategy. Options are
%                  'ordered' :uses them in the input order, but divided
%                  'random'  : orders them randomply
%                  'angularDistance': chooses the next subset with the
%                                     biggest angular distance with the ones used.
%
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
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------

%% Deal with input parameters

[blocksize,lambda,res,lamdbared,verbose,QualMeasOpts,OrderStrategy]=parse_inputs(proj,geo,angles,varargin);
measurequality=~isempty(QualMeasOpts);

if nargout>1
    computeL2=true;
else
    computeL2=false;
end

%% weigth matrices
% first order the projection angles
[alphablocks,orig_index]=order_subsets(angles,blocksize,OrderStrategy);


% Projection weigth, W
geoaux=geo;
geoaux.sVoxel(3)=max(geo.sDetector(2),geo.sVoxel(3)); % make sure lines are not cropped. One is for when image is bigger than detector and viceversa
geoaux.nVoxel=[2,2,2]'; % accurate enough?
geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'ray-voxel');  %
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;
% Back-Projection weigth, V
if ~isfield(geo,'mode')||~strcmp(geo.mode,'parallel')
    [x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
        -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
    A = permute(angles+pi/2, [1 3 2]);
    V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
    V=single(V);
else
    V=ones([geo.nVoxel(1:2).',length(angles)],'single');
end

clear A x y dx dz;



%% Iterate
errorL2=[];
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;



% TODO : Add options for Stopping criteria
for ii=1:niter
    
    % If verbose, time the algorithm
    if (ii==1 && verbose==1);tic;end
    % If quality is going to be measured, then we need to save previous image
    % THIS TAKES MEMORY!
    if measurequality
        res_prev=res;
    end
    
    
    for jj=1:length(alphablocks);
        % Get offsets
        if size(offOrigin,2)==length(angles)
            geo.offOrigin=offOrigin(:,orig_index{jj});
        end
        if size(offDetector,2)==length(angles)
            geo.offDetector=offDetector(:,orig_index{jj});
        end
        
        %proj is data: b=Ax
        %res= initial image is zero (default)
%         proj_err=proj(:,:,orig_index{jj})-Ax(res,geo,alphablocks{jj},'interpolated'); %                                 (b-Ax)
%         weighted_err=W(:,:,orig_index{jj}).*proj_err;                                 %                          W^-1 * (b-Ax)
%         backprj=Atb(weighted_err,geo,alphablocks{jj},'FDK');                          %                     At * W^-1 * (b-Ax)
%         weigth_backprj=bsxfun(@times,1./sum(V(:,:,orig_index{jj}),3),backprj);        %                 V * At * W^-1 * (b-Ax)
%         res=res+lambda*weigth_backprj;                                                % x= x + lambda * V * At * W^-1 * (b-Ax)
        
                res=res+lambda* bsxfun(@times,1./sum(V(:,:,orig_index{jj}),3),Atb(W(:,:,orig_index{jj}).*(proj(:,:,orig_index{jj})-Ax(res,geo,alphablocks{jj})),geo,alphablocks{jj}));
        
        
        % Non-negativity constrain
        res(res<0)=0;
        
        
    end
    
    % If quality is being measured
    if measurequality
        
        %Can save quality measure for every iteration here
        %See if some image quality measure should be used for every
        %iteration?
        qualMeasOut(:,ii)=Measure_Quality(res_prev,res,QualMeasOpts);
    end
    
    % reduce hyperparameter
    lambda=lambda*lamdbared;
    if computeL2
        % Compute error norm2 of b-Ax
        geo.offOrigin=offOrigin;
        geo.offDetector=offDetector;
        errornow=im3Dnorm(proj-Ax(res,geo,angles,'ray-voxel'),'L2');
        %     If the error is not minimized
        if ii~=1 && errornow>errorL2(end) % This 1.1 is for multigrid, we need to focus to only that case
            if verbose
                disp(['Convergence criteria met, exiting on iteration number:', num2str(ii)]);
            end
            return;
        end
        %     Store Error
        errorL2=[errorL2 errornow];
    end
    % If timing was asked
    if ii==1 && verbose==1
        expected_time=toc*(niter-1);
        expected_duration=toc*(niter);
        disp('OS-SART');
        disp(['Expected duration  :    ',secs2hms(expected_duration)]);
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
function [block_size,lambda,res,lamdbared,verbose,QualMeasOpts,OrderStrategy]=parse_inputs(proj,geo,alpha,argin)
opts=     {'BlockSize','lambda','Init','InitImg','Verbose','lambdaRed','QualMeas','OrderStrategy'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('CBCT:OS-SART:InvalidInput','Invalid number of inputs')
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
            if ~is2014bOrNewer
                warning('Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
            % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=0.99;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid lambda')
                end
                lambda=val;
            end
        case 'lambdaRed'
            if default
                lamdbared=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid lambda')
                end
                lamdbared=val;
            end
        case 'BlockSize'
            if default
                block_size=20;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid BlockSize')
                end
                block_size=val;
            end
            
        case 'Init'
            res=[];
            if default || strcmp(val,'none')
                res=zeros(geo.nVoxel','single');
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
                error('CBCT:OS_SART_CBCT:InvalidInput','Invalid Init option')
            end
            % % % % % % % ERROR
        case 'InitImg'
            if default
                continue;
            end
            if exist('initwithimage','var');
                if isequal(size(val),geo.nVoxel');
                    res=single(val);
                else
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid image for initialization');
                end
            end
        case 'QualMeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('CBCT:OS_SART_CBCT:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'OrderStrategy'
            if default
                OrderStrategy='angularDistance';
            else
                OrderStrategy=val;
            end
        otherwise
            error('CBCT:OS_SART_CBCT:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in OS_SART_CBCT()']);
    end
end

end




