function proj=addCTnoise(proj,varargin)
%ADDCTNOISE adds realistic noise to CT projections
%   addCTnoise(PROJ):  adds Poisson and Gaussian noise to the input data.
%
%   addCTnoise(PROJ,OPTS,VAL,..): adds noise with OPT options.
%                                 Possible options are :
%   'Poisson' : changes the average photon count for Poisson noise. Default
%               value is 1e5
%   'Gaussian': changes the mean and standard deviation of the electronic
%               Gaussian noise. Default value is [0 10]
%
% The computation of the noise projection follows
% Adaptive-weighted total variation minimization for sparse data toward low-dose x-ray computed tomography image reconstruction
% doi: 10.1088/0031-9155/57/23/7923
% and references in that work.
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
%% parse inputs


opts=     {'Poisson','Gaussian'};
defaults= [   1  ,  1 ];

% Check inputs
nVarargs = length(varargin);
if mod(nVarargs,2)
    error('CBCT:addnoise:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,varargin{ii}));
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
            ind=find(isequal(opt,varargin{jj}));
            jj=jj+1;
        end
        val=varargin{jj};
    end
    
    switch opt
        case 'Poisson'
            if default
                I0=1e5;
            else
               if ~isscalar(val);error('CBCT:addnoise:WorngInput','Input to Poisson should be scalar');end
               I0=val;
            end

        case 'Gaussian'
            if default
                m=0;
                sigma=10;
            else
                if (size(val,2)~=2);error('CBCT:addnoise:WorngInput','Input to Gaussian should be 1x2');end
                m=val(1);
                sigma=val(2);
            end
        otherwise
            error('CBCT:addnoise:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in OS_SART_CBCT()']);
    end
end

%% Add the noise
%//Lambert-Beer
Im=I0*exp(-proj);

% Photon noise + electronic noise
if areTheseToolboxesInstalled({'MATLAB','Statistics Toolbox'}) || areTheseToolboxesInstalled({'MATLAB','Statistics and Machine Learning Toolbox'})
    Im=poissrnd(Im)+randn(size(Im)).*sigma + m;
else
     warning(['You dont have Statistic toolbox, so poisson random noise is not available in MATLAB.',...
         java.lang.System.getProperty('line.separator').char,...
         'If you want to add that noise, use the following command:',...
         java.lang.System.getProperty('line.separator').char,...
         'Im=poissonrandom(I0*exp(-proj)); Im(Im<0)=1e-6; proj=single(log(I0./Im));',...
         java.lang.System.getProperty('line.separator').char,...,
         'With I0 ~ 10000'])
     Im=randn(size(Im)).*sigma + m; % this one is slower
end
Im(Im<0)=1e-6;
proj=single(log(I0./Im));
end