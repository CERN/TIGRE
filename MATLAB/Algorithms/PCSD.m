function [ f,qualMeasOut] = PCSD(proj,geo,angles,maxiter,varargin)
%PCSD solves the reconstruction problem using projection-controlled steepest descent method
%
%   PCSD(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem using
%   the projection data PROJ taken over ALPHA angles, corresponding to the
%   geometry described in GEO, using NITER iterations.
%
%   PCSD(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%
%   'lambda':      Sets the value of the hyperparameter for the SART iterations.
%                  Default is 1
%
%   'lambdared':   Reduction of lambda.Every iteration
%                  lambda=lambdared*lambda. Default is 0.99
%
%       'init':    Describes different initialization techniques.
%                   •  'none'     : Initializes the image to zeros (default)
%                   •  'FDK'      : initializes image to FDK reconstruction
%
%   'TViter':      Defines the amount of TV iterations performed per SART
%                  iteration. Default is 20
%
%   'maxL2err'     Maximum L2 error to accept an image as valid. This
%                  parameter is crucial for the algorithm, determines at
%                  what point an image should not be updated further.
%                  Default is 20% of the FDK L2 norm.
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
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
[beta,beta_red,f,ng,verbose,epsilon,QualMeasOpts,nonneg,gpuids,redundancy_weights,gt]=parse_inputs(proj,geo,angles,varargin);

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
qualMeasOut=zeros(length(QualMeasOpts),niter);

% does detector rotation exists?
if ~isfield(geo,'rotDetector')
    geo.rotDetector=[0;0;0];
end

%% Create weighting matrices for the SART step
% the reason we do this, instead of calling the SART function is not to
% recompute the weights every AwASD-POCS iteration, thus effectively doubling
% the computational time
% Projection weight, W

geoaux=geo;
geoaux.sVoxel([1 2])=geo.sVoxel([1 2])*1.1; % a Bit bigger, to avoid numerical division by zero (small number)
geoaux.sVoxel(3)=max(geo.sDetector(2),geo.sVoxel(3)); % make sure lines are not cropped. One is for when image is bigger than detector and viceversa
geoaux.nVoxel=[2,2,2]'; % accurate enough?
geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'Siddon','gpuids',gpuids);
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;

% Compute V
V=computeV(geo,angles,num2cell(angles),num2cell(1:length(angles)),'gpuids',gpuids);

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

%Initialize image.
%f=zeros(geo.nVoxel','single');

iter=0;
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;
rotDetector=geo.rotDetector;
stop_criteria=0;
DSD=geo.DSD;
DSO=geo.DSO;
%%
while ~stop_criteria %POCS
    % If quality is going to be measured, then we need to save previous image
    if measurequality && ~strcmp(QualMeasOpts,'error_norm')
        res_prev = f; % only store if necessary
    end
    if (iter==0 && verbose==1);tic;end
    iter=iter+1;
    
    %Estimation error in the projection domain
    est_proj=Ax(f,geo,angles,'interpolated','gpuids',gpuids);
    delta_p=im3Dnorm(est_proj-proj,'L2');
    
    %Enforcing ART along all projections if squared delta_p > epsilon
    if (delta_p^2)>epsilon
        for jj=1:size(angles,2)
            if size(offOrigin,2)==size(angles,2)
                geo.offOrigin=offOrigin(:,jj);
            end
            if size(offDetector,2)==size(angles,2)
                geo.offDetector=offDetector(:,jj);
            end
            if size(rotDetector,2)==size(angles,2)
                geo.rotDetector=rotDetector(:,jj);
            end
            if size(DSD,2)==size(angles,2)
                geo.DSD=DSD(jj);
            end
            if size(DSO,2)==size(angles,2)
                geo.DSO=DSO(jj);
            end
            f=f+beta* bsxfun(@times,1./V(:,:,jj),Atb(W(:,:,jj).*(proj(:,:,jj)-Ax(f,geo,angles(:,jj),'gpuids',gpuids)),geo,angles(:,jj),'gpuids',gpuids));
            
        end
    end
    
    %Non-negativity projection on all pixels
    if nonneg
        f=max(f,0);
    end
    
    geo.offDetector=offDetector;
    geo.offOrigin=offOrigin;
    geo.DSD=DSD;
    geo.DSO=DSO;
    geo.rotDetector=rotDetector;
    if measurequality
        qualMeasOut(:,iter)=Measure_Quality(res_prev,f,QualMeasOpts);
    end
    
    % Compute L2 error of actual image. Ax-b
    dd=im3Dnorm(Ax(f,geo,angles,'gpuids',gpuids)-proj,'L2');
    % Compute change in the image after last SART iteration
    dp_vec=(f-f0);
    
    if iter==1
        step=1;
    else
        step=delta_p/delta_p_first;
    end
    f0=f;
    %  TV MINIMIZATION
    % =========================================================================
    %  Call GPU to minimize TV
    f=minimizeTV(f0,step,ng,'gpuids',gpuids);    %   This is the MATLAB CODE, the functions are sill in the library, but CUDA is used nowadays
    %                                             for ii=1:ng
    %                                                 %delta=-0.00038 for thorax phantom
    %                                                 df=weighted_gradientTVnorm(f,delta);
    %                                                 df=df./im3Dnorm(df,'L2');
    %                                                 f=f-(step.*df);
    %                                             end
    
    % Compute change by TV min
    dg_vec=(f-f0);
    
    if iter==1
        delta_p_first=im3Dnorm((Ax(f0,geo,angles,'interpolated','gpuids',gpuids))-proj,'L2');
    end
    
    % Reduce SART step
    beta=beta*beta_red;
    
    % Check convergence criteria
    % ==========================================================================
    
    %Define c_alpha as in equation 21 in the journal
    c=dot(dg_vec(:),dp_vec(:))/(norm(dg_vec(:),2)*norm(dp_vec(:),2));
    %This c is examined to see if it is close to -1.0
    
    if (c<-0.99 && dd<=epsilon) || beta<0.005|| iter>maxiter
        if verbose
            disp(['Stopping criteria met']);
            disp(['   c    = ' num2str(c)]);
            disp(['   beta = ' num2str(beta)]);
            disp(['   iter = ' num2str(iter)]);
        end
        stop_criteria=true;
    end
    
    if (iter==1 && verbose==1)
        expected_time=toc*maxiter;
        disp('PCSD');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
end
end


function [beta,beta_red,f0,ng,verbose,epsilon,QualMeasOpts,nonneg,gpuids,redundancy_weights]=parse_inputs(proj,geo,angles,argin)
opts=     {'lambda','lambda_red','init','tviter','verbose','maxl2err','qualmeas','nonneg','gpuids','redundancy_weighting'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('CBCT:PCSD:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('CBCT:PCSD:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('CBCT:PCSD:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
                warning('Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
        % Lambda
        %  =========================================================================
        case 'lambda'
            if default
                beta=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:PCSD:InvalidInput','Invalid lambda')
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
                    error('TIGRE:PCSD:InvalidInput','Invalid lambda')
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
                    error('TIGRE:PCSD:InvalidInput','Invalid init')
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
        %  Maximum L2 error to have a "good image"
        %  =========================================================================
        case 'maxl2err'
            if default
                epsilon=im3Dnorm(FDK(proj,geo,angles))*0.2; %heuristic
            else
                epsilon=val;
            end
        %Image Quality Measure
        %  =========================================================================
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:PCSD:InvalidInput','Invalid quality measurement parameters');
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
        case 'redundancy_weighting'
            if default
                redundancy_weights = true;
            else
                redundancy_weights = val;
            end
        otherwise
            error('TIGRE:PCSD:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in PCSD()']);
            
    end
end

end
