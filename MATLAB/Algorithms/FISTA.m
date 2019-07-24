function [res] = FISTA(proj,geo,angles,niter,varargin)
% FISTA is a quadratically converging algorithm.

% 'hyper': This parameter should approximate the largest 
%          eigenvalue in the A matrix in the equation Ax-b and Atb. 
%          Empirical tests show that for, the headphantom object:
%               
%               geo.nVoxel = [64,64,64]'    ,      hyper (approx=) 2.e8
%               geo.nVoxel = [512,512,512]' ,      hyper (approx=) 2.e4
%          Default: 2.e8
% 'tviter':  Number of iterations of Im3ddenoise to use. Default: 20
% 'lambda':  Multiplier for the tvlambda used, which is proportional to 
%            L (hyper). Default: 0.1
% 'verbose': get feedback or not. Default: 1
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
[verbose,tviter,hyper,lambda]=parse_inputs(proj,geo,angles,varargin);
res = zeros(geo.nVoxel','single');
x_rec = res;
L = hyper;
bm = 1/L;
t = 1;
for ii = 1:niter
    if (ii==1);tic;end
    % gradient descent step
    res = res + bm * 2 * Atb(proj - Ax(res,geo,angles, 'ray-voxel'), geo, angles, 'matched');
    lambdaforTV = 2* bm* lambda;
    x_recold = x_rec;
    x_rec = im3DDenoise(res,'TV',tviter,1/lambdaforTV);  
    told = t;
    t = ( 1+sqrt(1+4*t*t) ) / 2;
    res= x_rec + (told-1)/t * (x_rec - x_recold);
    if (ii==1)&&(verbose==1);
        expected_time=toc*niter;
        disp('FISTA');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
end

end
%% Parse inputs
function [verbose,tviter,hyper,lambda]=parse_inputs(proj,geo,angles,argin)
opts=     {'lambda','tviter','verbose','hyper'};
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
            % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=0.1;
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
            else
                lambda=val;
            end
        case 'hyper'
            if default
                hyper = 2.e8;
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
            else 
                hyper = val;
            end
        case 'tviter'
            if default
                tviter = 20;
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
            else
                tviter = val;
            end
               
        otherwise
            error('TIGRE:OS_SART:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in OS_SART_CBCT()']);
    end
end
end
