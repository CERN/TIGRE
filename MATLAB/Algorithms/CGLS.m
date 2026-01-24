function [x,resL2,qualMeasOut]= CGLS(proj,geo,angles,niter,varargin)
% CGLS solves the CBCT problem using the conjugate gradient least
% squares
%
%  CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%  CGLS(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%
%  'Init'    Describes different initialization techniques.
%             * 'none'     : Initializes the image to zeros (default)
%             * 'FDK'      : initializes image to FDK reconstruction
%             * 'multigrid': Initializes image by solving the problem in
%                            small scale and increasing it when relative
%                            convergence is reached.
%             * 'image'    : Initialization using a user specified
%                            image. Not recommended unless you really
%                            know what you are doing.
%  'InitImg'    an image for the 'image' initialization. Avoid.
%
%  'groundTruth'  an image as ground truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
%  'restart'  true or false. By default the algorithm will restart when
%             loss of orthogonality is found. 
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

[verbose,x,QualMeasOpts,gpuids,gt,restart]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    x0=gt;
    clear gt
end
if nargout<3 && measurequality
    warning("Image metrics requested but none catched as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),niter);

resL2=zeros(1,niter);

% //doi: 10.1088/0031-9155/56/13/004
iter=0;
remember=0;
while iter<niter
    r=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
    p=Atb(r,geo,angles,'matched','gpuids',gpuids);
    gamma=norm(p(:),2)^2;
    for ii=iter:niter
        iter=iter+1;
        if measurequality && ~strcmp(QualMeasOpts,'error_norm')
            x0 = x; % only store if necessary
        end
        if (iter==1 && verbose);tic;end
        
        q=Ax(p,geo,angles,'Siddon','gpuids',gpuids);
        alpha=gamma/norm(q(:),2)^2;
        x=x+alpha*p;
        
        
        
        if measurequality
            qualMeasOut(:,iter)=Measure_Quality(x0,x,QualMeasOpts);
        end
        
        % The following should never happen, but the reality is that if we use
        % the residual from the algorithm, it starts diverging from this explicit residual value.
        % This is an interesting fact that I believe may be caused either by
        % the mismatch of the backprojection w.r.t the real adjoint, or
        % numerical issues related to doing several order of magnitude
        % difference operations on single precision numbers.
        aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
        resL2(iter)=im3Dnorm(aux,'L2');
        if iter>1 && resL2(iter)>resL2(iter-1)
            % we lost orthogonality, lets restart the algorithm unless the
            % user asked us not to. 
            
            % undo bad step. 
            x=x-alpha*p; 
            % if the restart didn't work. 
            if remember==iter || ~restart
                disp(['Algorithm stoped in iteration ', num2str(iter),' due to loss of orthogonality.'])
                return;
            end
            remember=iter;
            iter=iter-1;
            if verbose
               disp(['Orthogonality lost, restarting at iteration ', num2str(iter) ])
            end
            break  
            
        end
        % If step is adequate, then continue withg CGLS
        r=r-alpha*q;
        s=Atb(r,geo,angles,'matched','gpuids',gpuids);
        
        gamma1=norm(s(:),2)^2;
        beta=gamma1/gamma;
        gamma=gamma1;
        p=s+beta*p;
        
        
        if (iter==1 && verbose)
            expected_time=toc*niter;
            disp('CGLS');
            disp(['Expected duration   :    ',secs2hms(expected_time)]);
            disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
            disp('');
        end
    end
end
end


%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids,gt,restart]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','groundtruth','restart'};
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
            error('TIGRE:CGLS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end