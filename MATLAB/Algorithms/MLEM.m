function [res,qualMeasOut]=MLEM(proj,geo,angles,niter,varargin)
%MLEM solves the tomographic problem by using Maximum Likelihood Expectation
% Maximisation algorithm. 
%
%   MLEM(PROJ,GEO,ALPHA,NITER,opt) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry described in GEO, using NITER iterations.
%
% 'verbose': Get feedback or not. Default: 1
%
% 'init':    Describes different initialization techniques.
%             •  'none'     : Initializes the image to ones (default)
%             •  'FDK'      : Initializes image to FDK reconstruction
%
% 'QualMeas':    Asks the algorithm for a set of quality measurement
%                parameters. Input should contain a cell array of desired
%                quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                These will be computed in each iteration.
%
% 'groundTruth':  An image as ground truth, to be used if quality measures
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
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------
[verbose,res,QualMeasOpts,gpuids,gt]=parse_inputs(proj,geo,angles,varargin);
measurequality=~isempty(QualMeasOpts) | ~any(isnan(gt(:)));
if ~any(isnan(gt(:)))
    QualMeasOpts{end+1}='error_norm';
    res_prev=gt;
    clear gt
end
if nargout<2 && measurequality
    warning("Image metrics requested but none caught as output. Call the algorithm with 3 outputs to store them")
    measurequality=false;
end
qualMeasOut=zeros(length(QualMeasOpts),niter);


res = max(res,0);
% Back-projection weight, V
V = Atb(ones(size(proj),'single'),geo,angles,'matched','gpuids',gpuids);
V(V<=0.) = inf;

for ii=1:niter
    if measurequality && ~strcmp(QualMeasOpts,'error_norm')
        res_prev = res; % only store if necessary
    end
    if (ii==1);tic;end

    den = Ax(res,geo,angles,'gpuids',gpuids);
    den(den<=0.)=inf;
    
    imgupdate = Atb(proj./den, geo,angles,'matched','gpuids',gpuids)./V;
    res = max(res.*imgupdate,0.);
    
    if measurequality
        qualMeasOut(:,ii)=Measure_Quality(res_prev,res,QualMeasOpts);
    end
    
    if (ii==1)&&(verbose==1)
        expected_time=(toc)*niter;
        disp('MLEM');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end

end
end

%% Parse inputs
function [verbose,f0,QualMeasOpts,gpuids,gt]=parse_inputs(proj,geo,angles,argin)
opts = {'verbose','init','qualmeas','gpuids','groundtruth'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:MLEM:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:MLEM:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:MLEM:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
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
                f0=ones(geo.nVoxel','single');
            else
                if strcmp(val,'FDK')
                    f0=FDK(proj, geo, angles);
                else
                    error('TIGRE:MLEM:InvalidInput','Invalid init')
                end
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
                    error('TIGRE:MLEM:InvalidInput','Invalid quality measurement parameters');
                end
            end
        % GPUIDS
        % =========================================================================
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
        otherwise
            error('TIGRE:MLEM:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in MLEM()']);
    end
end
end

