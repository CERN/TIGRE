function [res,qualMeasOut] = FISTA(proj,geo,angles,niter,varargin)
% FISTA is a quadratically converging algorithm, modified from FISTA.
%
% It is based on the lazy-start FISTA modification in the following work:
% J. Liang and C. Schonlieb, "Faster FISTA", in 26th European Signal
% Processing Conference, 2018, pp.732-736. Lazy-start FISTA_mod
%

% 'hyper': This parameter should approximate the largest 
%          eigenvalue in the A matrix in the equation Ax-b and Atb. 
%          Empirical tests show that for, the headphantom object:
%               
%               geo.nVoxel = [64,64,64]'    ,      hyper (approx=) 2.e8
%               geo.nVoxel = [512,512,512]' ,      hyper (approx=) 2.e4
%          Default: 2.e8
% 'init':    Describes diferent initialization techniques.
%             •  'none'     : Initializes the image to zeros (default)
%             •  'FDK'      : intializes image to FDK reconstrucition
% 'tviter':  Number of iterations of Im3ddenoise to use. Default: 20
% 'lambda':  Multiplier for the tvlambda used, which is proportional to 
%            L (hyper). Default: 0.1
% 'fista_p': P parameter for "faster" FISTA (say 1/50). Default: 1 (standard
%            FISTA)
% 'fista_q': Q parameter for "faster" FISTA (say 1/10). Default: 1 (standard
%            FISTA)
% 'verbose': get feedback or not. Default: 1
%
% 'QualMeas'     Asks the algorithm for a set of quality measurement
%                parameters. Input should contain a cell array of desired
%                quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                These will be computed in each iteration.

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
% Coded by:           Ander Biguri, Reuben Lindroos
%--------------------------------------------------------------------------
[verbose,res,tviter,hyper,lambda,p,q,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,varargin);
%res = zeros(geo.nVoxel','single');
measurequality=~isempty(QualMeasOpts);

qualMeasOut=zeros(length(QualMeasOpts),niter);

x_rec = res;
L = hyper;
bm = 1/L;
t = 1;

for ii = 1:niter
    res0 = res;
    if (ii==1);tic;end
    % gradient descent step
    res = res + bm * 2 * Atb(proj - Ax(res,geo,angles, 'Siddon', 'gpuids', gpuids), geo, angles, 'matched', 'gpuids', gpuids);
    lambdaforTV = 2* bm* lambda;
    x_recold = x_rec;
    x_rec = im3DDenoise(res,'TV',tviter,1/lambdaforTV, 'gpuids', gpuids);  
    told = t;
    t = ( p+sqrt(q+r*t*t) ) / 2;
    res= x_rec + (told-1)/t * (x_rec - x_recold);

    if measurequality
        qualMeasOut(:,ii)=Measure_Quality(res0,res,QualMeasOpts);
    end

    if (ii==1)&&(verbose==1)
        expected_time=toc*niter;
        disp('FISTA');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
end

end
%% Parse inputs
function [verbose,f0,tviter,hyper,lambda,fista_p,fista_q,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,argin)
opts = {'lambda','init','tviter','verbose','hyper','fista_p','fista_q','qualmeas','gpuids'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:FISTA:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:FISTA:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:FISTA:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
        % Initial image
        %  =========================================================================
        case 'init'
            if default || strcmp(val,'none')
                f0=zeros(geo.nVoxel','single');
            else
                if strcmp(val,'FDK')
                    f0=FDK(proj, geo, angles);
                else
                    error('TIGRE:FISTA:InvalidInput','Invalid init')
                end
            end
        % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=0.1;
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:FISTA:InvalidInput','Invalid lambda')
            else
                lambda=val;
            end
        % hyperparameter
        % ==========================================================
        case 'hyper'
            if default
                hyper = 2.e8;
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:FISTA:InvalidInput','Invalid lambda')
            else 
                hyper = val;
            end
        % Number of iterations of TV
        %  =========================================================================
        case 'tviter'
            if default
                tviter = 20;
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:FISTA:InvalidInput','Invalid lambda')
            else
                tviter = val;
            end
        % FISTA parameter p
        %  =========================================================================
        case 'fista_p'
            if default
                fista_p = 1; % standard FISTA, 1/50 for "faster" FISTA
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:FISTA:InvalidInput','Invalid lambda')
            else
                fista_q = val;
            end
        % Number of iterations of TV
        %  =========================================================================
        case 'fista_q'
            if default
                fista_q = 1; % standard FISTA, 1/20 for "faster" FISTA
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:FISTA:InvalidInput','Invalid lambda')
            else
                fista_q = val;
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
                    error('TIGRE:FISTA:InvalidInput','Invalid quality measurement parameters');
                end
            end
        % GPU IDs
        %  =========================================================================
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
        otherwise
            error('TIGRE:FISTA:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in FISTA()']);
    end
end
end
