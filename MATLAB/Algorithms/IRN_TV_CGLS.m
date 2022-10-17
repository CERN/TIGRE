function [x_out,resL2,qualMeasOut]= IRN_TV_CGLS(proj,geo,angles,niter,varargin)
% IRN_TV_CGLS solves the IRN_TV_CGLS problem using the conjugate gradient least
% squares with Total Variation regularization, using an inner outer scheme
%
%  IRN_TV_CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%  IRN_TV_CGLS(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%  'lambda'  Value of regularization parameter lambda, default 10.
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
%  'groundTruth'  an image as grounf truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
%  'restart'  true or false. By default the algorithm will restart when
%             loss of ortogonality is found.
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

% % An Iteratively Reweighted Norm Algorithm for Total Variation Regularization
% % Paul Rodriguez and Brendt Wohlberg
% % 10.1109/ACSSC.2006.354879

[verbose,x0,QualMeasOpts,gpuids,lambda,niter_outer,gt,restart]=parse_inputs(proj,geo,angles,varargin);
x=x0;

niter_break=round(niter/niter_outer);

measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    x_prev=gt;
    clear gt
end
if nargout<3 && measurequality
    warning("Image metrics requested but none catched as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),niter);
resL2 = zeros(1,niter);

remember=[0];
iter=0;
% for iii = 1:niter_outer
while iter<niter
    if (iter==0 && verbose);tic;end
    
    % weights are N1 x N2 x N3
    W = build_weights(x);
    
    % Malena: We could avoid re-doing this, saving the initialisation
    % Just remember that we need cold restarts
    x=x0;
    
    prox_aux_1 = Ax(x,geo,angles,'Siddon','gpuids',gpuids);
    prox_aux_2 = Lx(W, x)*sqrt(lambda);
    
    r_aux_1 = proj - prox_aux_1;
    r_aux_2 = -prox_aux_2;
    
    p_aux_1 = Atb(r_aux_1 ,geo,angles,'matched','gpuids',gpuids);
    p_aux_2 = sqrt(lambda)*Ltx (W, r_aux_2);
    p = p_aux_1 + p_aux_2;
    
    gamma=norm(p(:),2)^2;
    
    while iter<niter %this will be broken by a parameter
        iter=iter+1;
        if measurequality && ~strcmp(QualMeasOpts,'error_norm')
            x_prev = x; % only store if necesary
        end
        q_aux_1 = Ax(p ,geo,angles,'Siddon','gpuids',gpuids);
        q_aux_2 = Lx (W, p)*sqrt(lambda);
        
        % q = cat(3, q_aux_1, q_aux_2{1},q_aux_2{2},q_aux_2{3}); % Probably never need to actually do this
        % alpha=gamma/norm(q(:),2)^2;
        alpha=gamma/(norm(q_aux_1(:))^2 + im3Dnorm(q_aux_2(:,:,:,1),'L2')^2 + im3Dnorm(q_aux_2(:,:,:,2),'L2')^2 + im3Dnorm(q_aux_2(:,:,:,3),'L2')^2);
        
        
        x=x+alpha*p;
        x_out=x;
        
        if measurequality
            qualMeasOut(:,iter)=Measure_Quality(x_prev,x,QualMeasOpts);
        end
        
        % The following should never happen, but the reallity is that if we use
        % the residual from the algorithm, it starts diverging from this explicit residual value.
        % This is an interesting fact that I believe may be caused either by
        % the mismatch of the backprojection w.r.t the real adjoint, or
        % numerical issues related to doing several order of magnitude
        % difference operations on single precission numbers.
        aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
        resL2(iter)=im3Dnorm(aux,'L2');
        if mod(iter,niter_break)~=1 && resL2(iter)>resL2(iter-1)
            % we lost orthogonality, lets restart the algorithm unless the
            % user asked us not to.
            iter
            niter_break
            % undo bad step.
            x=x-alpha*p;
            % if the restart didn't work.
            if remember==iter || ~restart
                disp(['Algorithm stoped in iteration ', num2str(iter),' due to loss of ortogonality.'])
                return;
            end
            remember=iter;
            iter=iter-1;
            if verbose
                disp(['Orthogonality lost, restarting at iteration ', num2str(iter) ])
            end
            break
        end
        
        % If step is adecuate, then continue withg CGLS
        r_aux_1 = r_aux_1-alpha*q_aux_1;
        r_aux_2=r_aux_2-alpha*q_aux_2;
        
        s_aux_1 = Atb(r_aux_1 ,geo,angles,'matched','gpuids',gpuids);
        s_aux_2 =  sqrt(lambda) * Ltx (W, r_aux_2);
        s = s_aux_1 + s_aux_2;
        
        gamma1=norm(s(:),2)^2;
        beta=gamma1/gamma;
        gamma=gamma1;
        p=s+beta*p;
        
        % replaces double explicit loops. Breaks the inner loop when its
        % time to do a cold restart.
        if mod(iter,niter_break)==0
            break
        end
        if (iter==1 && verbose)
            expected_time=toc*niter;
            disp('IRN_TV_CGLS');
            disp(['Expected duration   :    ',secs2hms(expected_time)]);
            disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
            disp('');
        end
    end
    
    
end
end
% % % Non-sense now, just giving output of the right dimensions
% % function out = Lx (sizeD3_1, sizeD3_2, x)
% % % L = speye(3*sizeD3_1,sizeD3_2 );
% % if prod(size(x)) == sizeD3_2 && 3*prod(size(x(:,:,1:end-1))) == sizeD3_1
% %     out = cat(3,x(:,:,1:end-1),x(:,:,1:end-1),x(:,:,1:end-1));
% % else
% %     error('wrong dimensions')
% % end
% % end
% %
% % function out = Ltx (sizeD3_1, sizeD3_2, x)
% % % L = speye(3*sizeD3_1,sizeD3_2 );
% % if prod(size(x)) == sizeD3_1
% %     out = x(:,:,1:512);
% % else
% %     error('wrong dimensions')
% % end
% % end

% How to make this efficient?

function W = build_weights ( x)
% Directional discrete derivatives
% Reflective boundary conditions
Dxx=zeros(size(x),'single');
Dyx=zeros(size(x),'single');
Dzx=zeros(size(x),'single');

Dxx(1:end-1,:,:)=x(1:end-1,:,:)-x(2:end,:,:);
Dyx(:,1:end-1,:)=x(:,1:end-1,:)-x(:,2:end,:);
Dzx(:,:,1:end-1)=x(:,:,1:end-1)-x(:,:,2:end);

W = (Dxx.^2+Dyx.^2+Dzx.^2+1e-6).^(-1/4); % Fix this...

end

function out = Lx (W, x)
% Directional discrete derivatives
% Reflective boundary conditions
Dxx=zeros(size(x),'single');
Dyx=zeros(size(x),'single');
Dzx=zeros(size(x),'single');

Dxx(1:end-1,:,:)=x(1:end-1,:,:)-x(2:end,:,:);
Dyx(:,1:end-1,:)=x(:,1:end-1,:)-x(:,2:end,:);
Dzx(:,:,1:end-1)=x(:,:,1:end-1)-x(:,:,2:end);
% Build weights - is it better to find the right rotation and add
% tensors?
out = cat(4,W .* Dxx, W .* Dyx, W .* Dzx);
end

function out = Ltx (W, x)
Wx_1 = W .* x(:,:,:,1);
Wx_2 = W .* x(:,:,:,2);
Wx_3 = W .* x(:,:,:,3);

% Left here, but this is how to make a transpose
DxtWx_1=Wx_1;
DytWx_2=Wx_2;
DztWx_3=Wx_3;

DxtWx_1(2:end-1,:,:)=Wx_1(2:end-1,:,:)-Wx_1(1:end-2,:,:);
DxtWx_1(end,:,:)=-Wx_1(end-1,:,:);

DytWx_2(:,2:end-1,:)=Wx_2(:,2:end-1,:)-Wx_2(:,1:end-2,:);
DytWx_2(:,end,:)=-Wx_2(:,end-1,:);

DztWx_3(:,:,2:end-1)=Wx_3(:,:,2:end-1)-Wx_3(:,:,1:end-2);
DztWx_3(:,:,end)=-Wx_3(:,:,end-1);

out = DxtWx_1 + DytWx_2 + DztWx_3;
end
%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids,lambda,niter_outer,gt,restart]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','lambda','niter_outer','groundtruth','restart'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:IRN_TV_CGLS:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:IRN_TV_CGLS:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:IRN_TV_CGLS:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
                error('TIGRE:IRN_TV_CGLS:InvalidInput','Invalid Init option')
            end
        case 'niter_outer'
            if default
                niter_outer=5;
            else
                niter_outer=val;
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
                    error('TIGRE:IRN_TV_CGLS:InvalidInput','Invalid image for initialization');
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
                    error('TIGRE:IRN_TV_CGLS:InvalidInput','Invalid quality measurement parameters');
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
                lambda=10;
            else
                lambda=val;
            end
        case 'groundtruth'
            if default
                gt=nan;
            else
                gt=val;
            end
        case 'restart'
            if default
                restart=true;
            else
                restart=val;
            end
        otherwise
            error('TIGRE:IRN_TV_CGLS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in IRN_TV_CGLS()']);
    end
end


end