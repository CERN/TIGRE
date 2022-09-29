function [x_out,errorL2,qualMeasOut]= IRN_TV_CGLS(proj,geo,angles,niter,varargin)
% CGLS solves the CBCT problem using the conjugate gradient least
% squares
% 
%  CGLS(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
% 
%  CGLS(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
% 
% 
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
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------

%%

% % An Iteratively Reweighted Norm Algorithm for Total Variation Regularization
% % Paul Rodriguez and Brendt Wohlberg
% % 10.1109/ACSSC.2006.354879
% 
% errorL2 = [];
% for i = 1:10 % Needs to change
% % Update weight
% Build L = W * D as an efficient handle. Think about this.
%
% % Call CGLS (Maybe try with LSQR just for fun)
% [x ,errorL2_i,qualMeasOut] = CGLS(proj,geo,angles,niter,varargin);
%
% errorL2 = [errorL2, errorL2_i];
% % Malena: not sure what qualMeasOut is, so I do not know if it should be stored or
% % updated after every inner cycle
% end

% % Start coding with an L just to test dimensions
%
% sizeD3_1 =    (geo.nVoxel(1)-1)*geo.nVoxel(1)*geo.nVoxel(1) ...
%             + geo.nVoxel(1)*(geo.nVoxel(1)-1)*geo.nVoxel(1) ...
%             + geo.nVoxel(1)*geo.nVoxel(1)*(geo.nVoxel(1)-1);
% sizeD3_2 =  geo.nVoxel(1)*geo.nVoxel(1)*geo.nVoxel(1);

lambda = 10; %% ||Ax-b|| + lambda R(x)

[verbose,x,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,varargin);

niter_outer = 15;

measurequality=~isempty(QualMeasOpts);
qualMeasOut=zeros(length(QualMeasOpts),niter*niter_outer);
errorL2=zeros(1,niter*niter_outer);

for iii = 1:niter_outer

    % weights are N1 x N2 x N3
    W = build_weights (x);

    % Malena: We could avoid re-doing this, saving the initialisation
    % Just remember that we need cold restarts
    [verbose,x,QualMeasOpts,gpuids] = parse_inputs(proj,geo,angles,varargin);

    prox_aux_1 = Ax(x,geo,angles,'Siddon','gpuids',gpuids);
    prox_aux_2 = Lx(W, x);
    prox_aux_2 = cellfun(@(x) x*(sqrt(lambda)),prox_aux_2,'un',0);

    r_aux_1 = proj - prox_aux_1;
    r_aux_2 = cellfun(@(x) x*(-1),prox_aux_2,'un',0);

    % Malena: changed the format, r_aux_2 is 3
    % r = cat(3,r_aux_1, r_aux_2); % Malena: size guide, erase later, N x N x (100 + N-1)

    p_aux_1 = Atb(r_aux_1 ,geo,angles,'matched','gpuids',gpuids);
    p_aux_2 = sqrt(lambda)*Ltx (W, r_aux_2);
    p = p_aux_1 + p_aux_2;

    gamma=norm(p(:),2)^2;

    for ii=1:niter
        x0 = x;
        if (ii==1 && verbose);tic;end

        q_aux_1 = Ax(p ,geo,angles,'Siddon','gpuids',gpuids);
        q_aux_2 = Lx (W, p);
        q_aux_2 = cellfun(@(x) x*sqrt(lambda),q_aux_2,'un',0);

        % q = cat(3, q_aux_1, q_aux_2{1},q_aux_2{2},q_aux_2{3}); % Probably never need to actually do this
        % alpha=gamma/norm(q(:),2)^2;
        alpha=gamma/(norm(q_aux_1(:))^2 + norm(q_aux_2{1}(:))^2 + norm(q_aux_2{2}(:))^2 + norm(q_aux_2{3}(:))^2);


        x=x+alpha*p;
        x_out{iii}=x;
        aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);   % expensive, is there any way to check this better?
        % residual norm or the original least squares (not Tikhonov).
        % Think if that is what we want of the NE residual

        errorL2((iii-1)*niter+ii)=im3Dnorm(aux,'L2');

        if measurequality
            qualMeasOut(:,(iii-1)*niter+ii)=Measure_Quality(x0,x,QualMeasOpts);
        end

        % This might happen now, we are minimizing the normal eqns residual
%         if ii>1 && errorL2(ii)>errorL2(ii-1)
%             % OUT!
%             x=x-alpha*p;
%             if verbose
%                 disp(['CGLS stoped in iteration N', num2str(ii),' due to divergence.'])
%             end
%             return;
%         end

        % If step is adecuate, then continue withg CGLS
        r_aux_1 = r_aux_1-alpha*q_aux_1;
        r_aux_2 = cellfun(@(a,b)a-alpha*b, r_aux_2, q_aux_2, 'uni',0);

        s_aux_1 = Atb(r_aux_1 ,geo,angles,'matched','gpuids',gpuids);
        s_aux_2 =  sqrt(lambda) * Ltx (W, r_aux_2);
        s = s_aux_1 + s_aux_2;

        gamma1=norm(s(:),2)^2;
        beta=gamma1/gamma;
        gamma=gamma1;
        p=s+beta*p;

        if (iii==1 && ii ==1 && verbose)
            expected_time=toc*niter*niter_outer;
            disp('CGLS');
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
    Dxx=cat(1,x(1:end-1,:,:)-x(2:end,:,:),zeros(1,size(x,2),size(x,3)));
    Dyx=cat(2,x(:,1:end-1,:)-x(:,2:end,:),zeros(size(x,1),1,size(x,3)));
    Dzx=cat(3,x(:,:,1:end-1)-x(:,:,2:end),zeros(size(x,1),size(x,2),1));

    W = (Dxx.^2+Dyx.^2+Dzx.^2+1e-6).^(-1/4); % Fix this...

end

function out = Lx (W, x)
    % Directional discrete derivatives
    % Reflective boundary conditions 
    Dxx=cat(1,x(1:end-1,:,:)-x(2:end,:,:),zeros(1,size(x,2),size(x,3)));
    Dyx=cat(2,x(:,1:end-1,:)-x(:,2:end,:),zeros(size(x,1),1,size(x,3)));
    Dzx=cat(3,x(:,:,1:end-1)-x(:,:,2:end),zeros(size(x,1),size(x,2),1));

    % Build weights - is it better to find the right rotation and add
    % tensors?
    WDxx = W .* Dxx;
    WDyx = W .* Dyx;
    WDzx = W .* Dzx;

    out = {WDxx WDyx WDzx};
end

function out = Ltx (W, x) % Input is a cell 
    Wx_1 = W .* x{1};
    Wx_2 = W .* x{2};
    Wx_3 = W .* x{3};

    % Left here, but this is how to make a transpose
    DxtWx_1=cat(1,Wx_1(1,:,:),Wx_1(2:end-1,:,:)-Wx_1(1:end-2,:,:),-Wx_1(end-1,:,:));
    DytWx_2=cat(2,Wx_2(:,1,:),Wx_2(:,2:end-1,:)-Wx_2(:,1:end-2,:),-Wx_2(:,end-1,:));
    DztWx_3=cat(3,Wx_3(:,:,1),Wx_3(:,:,2:end-1)-Wx_3(:,:,1:end-2),-Wx_3(:,:,end-1));

    out = DxtWx_1 + DytWx_2 + DztWx_3;
end
%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids'};
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
        otherwise 
            error('TIGRE:CGLS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end