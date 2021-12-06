function proj=addCTnoise(proj, varargin)
%ADDCTNOISE adds realistic noise to CT projections
%   addCTnoise(PROJ):  adds Poisson and Gaussian noise to the input data.
%
%   addCTnoise(PROJ,OPTS,VAL,..): adds noise with OPT options.
%                                 Possible options are :
%   'Poisson' : changes the average photon count for Poisson noise. Default
%               value is 1e5
%   'Gaussian': changes the mean and standard deviation of the electronic
%               Gaussian noise. Default value is [0 0.5]
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
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
%% parse inputs


opts=     {'Poisson','Gaussian', 'Implementation'};
defaults= [   1  ,  1  ,  1 ];

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
                I0=60000;
                if max(proj(:)>I0)
                    I0=max(proj(:))/5;
                end
            else
                if ~isscalar(val);error('CBCT:addnoise:WrongInput','Input to Poisson should be scalar');end
                I0=val;
            end
            
            
        case 'Gaussian'
            if default
                m=0;
                sigma=0.5;
            else
                if (size(val,2)~=2);error('CBCT:addnoise:WrongInput','Input to Gaussian should be 1x2');end
                m=val(1);
                sigma=val(2);
            end
        
        case 'Implementation'
            if default
                implementation='cuda';
            else
                implementation=val;
            end
            
        otherwise
            error('CBCT:addnoise:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in addCTnoise()']);
    end
end

%% Add the noise
%//Lambert-Beer
Im=I0*exp(-proj/max(proj(:)));

% Photon noise + electronic noise
if strcmp(implementation, 'matlab')
    if areTheseToolboxesInstalled({'MATLAB','Image Processing Toolbox'})
        %disp('Using Image Processing Toolbox');
        Im=imnoise(Im/I0,'poisson')*I0;
        Im=imnoise(Im/I0,'gaussian',m,sigma/I0)*I0;
    elseif areTheseToolboxesInstalled({'MATLAB','Statistics Toolbox'}) || areTheseToolboxesInstalled({'MATLAB','Statistics and Machine Learning Toolbox'})
        %disp('Using Statistics Toolbox');
        Im=poissrnd(Im)+randn(size(Im)).*sigma + m;
    else
        warning(['You don''t have Image Processing Toolbox nor Statistic Toolbox.',...
            java.lang.System.getProperty('line.separator').char,...
            'CUDA version is used.']);
        %disp('Using CUDA ..');
        Im=AddNoise(Im, m, sigma);
    end
else % if implementation == 'cuda'
    %disp('Using CUDA');
    Im=AddNoise(Im, m, sigma);
end
Im(Im<=0)=1e-6;
proj=single(log(I0./Im))*max(proj(:));
end