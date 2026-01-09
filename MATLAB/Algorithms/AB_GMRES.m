function [x,resL2,qualMeasOut]= AB_GMRES(proj,geo,angles,niter,varargin)

% AB_GMRES solves the CBCT problem using AB_GMRES.
% This is a stable Krylov method for when the backprojector is not adjoint
%
%  AB_GMRES(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
%  AB_GMRES(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
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
%  'groundTruth'  an image as ground truth, to be used if quality measures
%                 are requested, to plot their change w.r.t. this known
%                 data.
% 'backprojector' Describes which backprojector to use
%
%                *'FDK': This will use the FDK algorithm as a backprojector
%                        NOTE: not the backprojector TIGRE calls "FDK", but
%                        the actual algorithm
%                *'matched': (DEFAULT) uses the pseudo-matched backprojectors 
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
% Coded by:           Malena Sabate Landman, Ander Biguri
%--------------------------------------------------------------------------

%%

[verbose,x,QualMeasOpts,gpuids,gt,bp]=parse_inputs(proj,geo,angles,varargin);

if strcmpi(bp,'FDK')
    backproject=@(proj,geo,angles,gpuids)FDK(proj,geo,angles,'gpuids',gpuids);
else
    backproject=@(proj,geo,angles,gpuids)Atb(proj,geo,angles,'matched','gpuids',gpuids);
end

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
if nargout==1
    warning("Divergence caused by numerical issues not checked, due to its high computational cost for AB-GMRES")
else
    warning("AB-GMRES takes very long to compute residual norms and quality measures and error norms");
end

qualMeasOut=zeros(length(QualMeasOpts),niter);
resL2=zeros(1,niter);


% Per Cristian Hansen:
% GMRES methods for tomographic reconstruction with an unmatched back projector

w=zeros(numel(proj),niter+1,'single');
r=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
w(:,1) = r(:)/norm(r(:),2);
h=zeros(niter+1,niter);
for k=1:niter
    if measurequality && ~strcmp(QualMeasOpts,'error_norm')
        x0 = x; % only store if necessary
    end
    if (k==1 && verbose);tic;end
    
    qk=Ax(backproject(reshape(w(:,k),geo.nDetector(1),geo.nDetector(2),length(angles)),geo,angles,gpuids),geo,angles,'Siddon','gpuids',gpuids);
    e1=zeros(k+1,1);
    e1(1)=1;
    for ii=1:k
        h(ii,k)=sum(qk(:).*w(:,ii));
        qk(:)=qk(:)-h(ii,k)*w(:,ii);
    end
    h(k+1,k)=norm(qk(:),2);
    w(:,k+1)=qk(:)/h(k+1,k);
    y=h(1:k+1,1:k)\(e1*norm(r(:),2));
    if measurequality
        qualMeasOut(:,k)=Measure_Quality(x0,compute_res(x,w(:,1:k),y,geo,angles,backproject,gpuids) ,QualMeasOpts);
    end

    if nargout>1
        resL2(k)=im3Dnorm(proj-Ax(compute_res(x,w(:,1:k),y,geo,angles,backproject,gpuids),geo,angles,'Siddon','gpuids',gpuids),'L2');
        if k>1 && resL2(k)>resL2(k-1)
            x=compute_res(x,w(:,1:end-1),y(1),geo,angles,backproject,gpuids);
            disp(['Algorithm stoped in iteration ', num2str(k),' due to loss of orthogonality.'])
            return;
        end
        
    end
    
    if (k==1 && verbose)
        expected_time=toc*niter;
        disp('AB-GMRES');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
end
x=compute_res(x,w(:,1:end-1),y,geo,angles,backproject,gpuids);



end

% x is not explicit in this algorith, and its quite expensive to compute
function x=compute_res(x,w,y,geo,angles,backproject,gpuIds)

for ii=1:size(w,2)
    x=x+backproject(reshape(w(:,ii),geo.nDetector(1),geo.nDetector(2),length(angles)),geo,angles,gpuIds)*y(ii);
end

end

%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids,gt,bp]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids','groundtruth','backprojector'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:AB-GMRES:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:AB-GMRES:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:AB-GMRES:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
                error('TIGRE:AB-GMRES:InvalidInput','Invalid Init option')
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
                    error('TIGRE:AB-GMRES:InvalidInput','Invalid image for initialization');
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
                    error('TIGRE:AB-GMRES:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE:AB-GMRES:Verbose mode not available for older versions than MATLAB R2014b');
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
        case 'backprojector'
            if default
                bp='matched';
            else
                bp=val;
            end
        otherwise
            error('TIGRE:AB-GMRES:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end