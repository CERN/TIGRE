function [ f,qualMeasOut ] = ASD_POCS(proj,geo,angles,maxiter,varargin)
%ASD_POCS Solves the ASD_POCS total variation constrained image in 3D
% tomography.
%
%   ASD_POCS(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%   ASD_POCS(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%
%   'lambda':      Sets the value of the hyperparameter for the SART iterations.
%                  Default is 1
%
%   'lambdared':   Reduction of lambda. Every iteration
%                  lambda = lambdared * lambda. Default is 0.99
%
%   'init':        Describes different initialization techniques.
%                   •  'none'     : Initializes the image to zeros (default)
%
%                   •  'FDK'      : Initializes image to FDK reconstruction
%   'TViter':      Defines the amount of TV iterations performed per SART
%                  iteration. Default is 20
%
%   'alpha':       Defines the TV hyperparameter. Default is 0.002
%
%   'alpha_red':   Defines the reduction rate of the TV hyperparameter
%
%   'Ratio':       The maximum allowed image/TV update ratio. If the TV
%                  update changes the image more than this, the parameter
%                  will be reduced. Default is 0.95
%   'maxL2err'     Maximum L2 error to accept an image as valid. This
%                  parameter is crucial for the algorithm, determines at
%                  what point an image should not be updated further.
%                  Default is 20% of the FDK L2 norm.
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%
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
% Coded by:           Ander Biguri and Manasavee Lohvithee
%--------------------------------------------------------------------------



%% parse inputs
blocksize=1;
[beta,beta_red,f,ng,verbose,alpha,alpha_red,rmax,epsilon,OrderStrategy,QualMeasOpts,nonneg,gpuids,redundancy_weights,gt]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    res_prev=gt;
    clear gt
end
if nargout<2 && measurequality
    warning("Image metrics requested but none catched as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),maxiter);


[alphablocks,orig_index]=order_subsets(angles,blocksize,OrderStrategy);
angles_reorder=cell2mat(alphablocks);
index_angles=cell2mat(orig_index);

% does detector rotation exists?
if ~isfield(geo,'rotDetector')
    geo.rotDetector=[0;0;0];
end
%% Create weighting matrices for the SART step
% the reason we do this, instead of calling the SART function is not to
% recompute the weights every ASD-POCS iteration, thus effectively doubling
% the computational time
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


%%
stop_criteria=0;
iter=0;
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;
rotDetector=geo.rotDetector;
DSD=geo.DSD;
DSO=geo.DSO;
while ~stop_criteria %POCS
    % If quality is going to be measured, then we need to save previous image
    if measurequality && ~strcmp(QualMeasOpts,'error_norm')
        res_prev = f; % only store if necessary
    end
    
    f0=f;
    if (iter==0 && verbose==1);tic;end
    iter=iter+1;
    
    for jj=1:size(angles,2)
        if size(offOrigin,2)==size(angles,2)
            geo.offOrigin=offOrigin(:,index_angles(:,jj));
        end
        if size(offDetector,2)==size(angles,2)
            geo.offDetector=offDetector(:,index_angles(:,jj));
        end
        if size(rotDetector,2)==size(angles,2)
            geo.rotDetector=rotDetector(:,index_angles(:,jj));
        end
        if size(DSD,2)==size(angles,2)
            geo.DSD=DSD(jj);
        end
        if size(DSO,2)==size(angles,2)
            geo.DSO=DSO(jj);
        end
        %         proj_err=proj(:,:,jj)-Ax(f,geo,angles(:,jj));          %                                 (b-Ax)
        %         weighted_err=W(:,:,jj).*proj_err;                   %                          W^-1 * (b-Ax)
        %         backprj=Atb(weighted_err,geo,angles(:,jj));            %                     At * W^-1 * (b-Ax)
        %         weight_backprj=bsxfun(@times,1./V(:,:,jj),backprj); %                 V * At * W^-1 * (b-Ax)
        %         f=f+beta*weight_backprj;                          % x= x + lambda * V * At * W^-1 * (b-Ax)
        % Enforce positivity
        f=f+beta* bsxfun(@times,1./V(:,:,jj),Atb(W(:,:,index_angles(:,jj)).*(proj(:,:,index_angles(:,jj))-Ax(f,geo,angles_reorder(:,jj),'gpuids',gpuids)),geo,angles_reorder(:,jj),'gpuids',gpuids));
        % non-negativity constrain
        if nonneg
            f=max(f,0);
        end
    end
    
    geo.offDetector=offDetector;
    geo.offOrigin=offOrigin;
    geo.DSD=DSD;
    geo.DSO=DSO;
    geo.rotDetector=rotDetector;
    % Save copy of image.
    if measurequality
        qualMeasOut(:,iter)=Measure_Quality(res_prev,f,QualMeasOpts);
    end
    % compute L2 error of actual image. Ax-b
    dd=im3Dnorm(Ax(f,geo,angles,'gpuids',gpuids)-proj,'L2');
    % compute change in the image after last SART iteration
    dp_vec=(f-f0);
    dp=im3Dnorm(dp_vec,'L2');
    
    if iter==1
        dtvg=alpha*dp;
        %Convert the steepest-descent step-size from a fraction of a
        %step-size to an absolute image distance on the first iteration.
    end
    f0=f;
    
    %  TV MINIMIZATION
    % =========================================================================
    %  Call GPU to minimize TV
    f=minimizeTV(f0,dtvg,ng,'gpuids',gpuids);    %   This is the MATLAB CODE, the functions are sill in the library, but CUDA is used nowadays
    %                                   for ii=1:ng
    % %                                 Steepest descend of TV norm
    %                                      tv(ng*(iter-1)+ii)=im3Dnorm(f,'TV','forward');
    %                                      df=weighted_gradientTVnorm2(f,0.002);
    %                                      df=df./im3Dnorm(df,'L2');
    %                                      f=f-dtvg.*df;
    %                                    end
    
    % update parameters
    % ==========================================================================
    
    % compute change by TV min
    dg_vec=(f-f0);
    dg=im3Dnorm(dg_vec,'L2');
    % if change in TV is bigger than the change in SART AND image error is still bigger than acceptable
    if dg>rmax*dp && dd>epsilon
        dtvg=dtvg*alpha_red;
    end
    % reduce SART step
    beta=beta*beta_red;
    % Check convergence criteria
    % ==========================================================================
    
    %Define c_alpha as in equation 21 in the journal
    c=dot(dg_vec(:),dp_vec(:))/(norm(dg_vec(:),2)*norm(dp_vec(:),2));
    %This c is examined to see if it is close to -1.0
    
    if (c<-0.99 && dd<=epsilon) || beta<0.005|| iter>=maxiter
        if verbose
            disp(['Stopping criteria met']);
            disp(['   c    = ' num2str(c),      '(Desired: c<-0.99)']);
            disp(['   beta = ' num2str(beta),   '(Desired: beta<0.005)']);
            disp(['   iter = ' num2str(iter), ]);
        end
        stop_criteria=true;
    end
    if (iter==1 && verbose==1)
        expected_time=toc*maxiter;
        disp('ADS_POCS');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
end



end

function [beta,beta_red,f0,ng,verbose,alpha,alpha_red,rmax,epsilon,OrderStrategy,QualMeasOpts,nonneg,gpuids,redundancy_weights,gt]=parse_inputs(proj,geo,angles,argin)

opts=     {'lambda','lambda_red','init','tviter','verbose','alpha','alpha_red','ratio','maxl2err','orderstrategy','qualmeas','nonneg','gpuids','redundancy_weighting','groundtruth'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:ASD_POCS:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:ASD_POCS:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extract value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end
        if isempty(ind)
            error('TIGRE:ASD_POCS:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
        end
        val=argin{jj};
    end
    % parse inputs
    switch opt
        % Verbose
        %  =========================================================================
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
            % Lambda
            %  =========================================================================
            % Its called beta in ASD-POCS
        case 'lambda'
            if default
                beta=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:ASD_POCS:InvalidInput','Invalid lambda')
                end
                beta=val;
            end
            % Lambda reduction
            %  =========================================================================
        case 'lambda_red'
            if default
                beta_red=0.99;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:ASD_POCS:InvalidInput','Invalid lambda')
                end
                beta_red=val;
            end
            % Initial image
            %  =========================================================================
        case 'init'
            if default || strcmp(val,'none')
                f0=zeros(geo.nVoxel','single');
                
            else
                if strcmp(val,'FDK')
                    f0=FDK(proj, geo, angles);
                else
                    error('TIGRE:ASD_POCS:InvalidInput','Invalid init')
                    
                end
            end
            % Number of iterations of TV
            %  =========================================================================
        case 'tviter'
            if default
                ng=20;
            else
                ng=val;
            end
            %  TV hyperparameter
            %  =========================================================================
        case 'alpha'
            if default
                alpha=0.002; % 0.2
            else
                alpha=val;
            end
            %  TV hyperparameter redution
            %  =========================================================================
        case 'alpha_red'
            if default
                alpha_red=0.95;
            else
                alpha_red=val;
            end
            %  Maximum update ratio
            %  =========================================================================
        case 'ratio'
            if default
                rmax=0.95;
            else
                rmax=val;
            end
            %  Maximum L2 error to have a "good image"
            %  =========================================================================
        case 'maxl2err'
            if default
                epsilon=im3Dnorm(Ax(FDK(proj,geo,angles),geo,angles)-proj,'L2')*0.2; %heuristic
            else
                epsilon=val;
            end
            %  Order strategy
            %  =========================================================================
        case 'orderstrategy'
            if default
                OrderStrategy='random';
            else
                OrderStrategy=val;
            end
            % Image Quality Measure
            %  =========================================================================
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:ASD_POCS:InvalidInput','Invalid quality measurement parameters');
                end
            end
            %  Non negative
            %  =========================================================================
        case 'nonneg'
            if default
                nonneg=true;
            else
                nonneg=val;
            end
            %  GPU Ids
            %  =========================================================================
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
            % Data Redundancy weighting
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
            error('TIGRE:ASD_POCS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in ASD_POCS()']);
            
    end
end

end


