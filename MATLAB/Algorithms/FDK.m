function [res]=FDK(proj,geo,angles,varargin)
%FDK solves Cone Beam CT image reconstruction using Feldkam Davis Kress
% algorithm (filtered backprojection) 
%
%   FDK(PROJ,GEO,ANGLES) solves the reconstruction problem
%   using the projection data PROJ taken over ANGLES angles, corresponding
%   to the geometry described in GEO.
%
%   FDK(PROJ,GEO,ANGLES,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%   'parker': adds parker weights for limited angle scans. Default TRUE
%
%   'wang': adds detector offset weights. Default TRUE
%
%   'filter': selection of filter. Default 'ram-lak' (ramp)
%              options are: 
%                  'ram-lak' (ramp)
%                  'shepp-logan'
%                  'cosine'
%                  'hamming'  
%                  'hann'
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
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Kyungsang Kim, modified by Ander Biguri, Brandon Nelson 
%--------------------------------------------------------------------------
[filter,parker,dowang,gpuids]=parse_inputs(proj,geo,angles,varargin);

geo=checkGeo(geo,angles);
geo.filter=filter;

if dowang
    disp('FDK: applying detector offset weights')
    % Zero-padding to avoid FFT-induced aliasing
    [zproj, zgeo, theta] = zeropadding(proj, geo);
    % Preweighting using Wang function to save memory
    [proj, ~] = preweighting2(zproj, zgeo, theta);
    
    %% Replace original proj and geo
    % proj = proj_w;
    geo = zgeo;
end

if size(geo.offDetector,2)==1
    offset=repmat(geo.offDetector,[1 size(angles,2)]);
else
    offset=geo.offDetector;
end

%% Weight
proj=permute(proj,[2 1 3]);
for ii=1:size(angles,2)
    
    us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1,ii);
    vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2,ii);
    [uu,vv] = meshgrid(us,vs); % detector
    
    % Create weight according to each detector element
    w = (geo.DSD(ii))./sqrt((geo.DSD(ii))^2+uu.^2 + vv.^2);
    
    % Multiply the weights with projection data
    proj(:,:,ii) = proj(:,:,ii).*w';
end
%% Fourier transform based filtering
proj = filtering(proj,geo,angles,parker); % Not sure if offsets are good in here

% RMFIELD Remove fields from a structure array.
geo=rmfield(geo,'filter');
%% backproject
%%%
% [proj, w] = preweighting(proj,geo);
% imshow(w,[])
%%%
res=Atb((proj),geo,angles, 'gpuids', gpuids); % Weighting is inside


end

function [zproj, zgeo, theta] = zeropadding(proj, geo)
% ZEROPADDING as preprocessing for preweighting
zgeo = geo;

padwidth = fix(2*geo.offDetector(1)./geo.dDetector(1));
zgeo.offDetector(1,:) = geo.offDetector(1,:) - padwidth/2 * geo.dDetector(1);
zgeo.nDetector(1) = abs(padwidth) + geo.nDetector(1);
zgeo.sDetector(1) = zgeo.nDetector(1) * zgeo.dDetector(1);

theta = (geo.sDetector(1)/2 - abs(geo.offDetector(1)))...
        * sign(geo.offDetector(1));
% Pad on the left size when offset >0
if(geo.offDetector(1)>0)
    for ii = 1:size(proj,3)
        zproj(:,:,ii) = [zeros(size(proj,1), padwidth), proj(:,:,ii)];
    end
else
    for ii = 1:size(proj,3)
        zproj(:,:,ii) = [proj(:,:,ii), zeros(size(proj,1), abs(padwidth))];
    end
end

end

function [proj_w, w] = preweighting(proj,geo,theta)
% Preweighting using Wang function
% Ref: 
%    Wang, Ge. X-ray micro-CT with a displaced detector array. Medical Physics, 2002,29(7):1634-1636.
offset = geo.offDetector(1);
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + abs(offset);

abstheta=abs(theta);

w = ones(size(proj(:,:,1)));
for ii = 1:geo.nDetector
    t = us(ii);
    if(abs(t) <= abstheta)
        w(:,ii) = 0.5*(sin((pi/2)*atan(t/geo.DSO(1))/(atan(theta/geo.DSO(1)))) + 1);
    end
    if(t<-abstheta)
        w(:,ii) = 0;
    end
end

if(theta<0)
    w = fliplr(w);
end

for ii = 1:size(proj,3)
    proj_w(:,:,ii) = proj(:,:,ii).*w*2;
end

end
function [proj_w, w] = preweighting2(proj,geo,theta)
% Preweighting using Wang function
% Ref: 
%    Wang, Ge. X-ray micro-CT with a displaced detector array. Medical Physics, 2002,29(7):1634-1636.
offset = geo.offDetector(1);
offset = offset + (geo.DSD(1)/geo.DSO(1))*geo.COR(1);
us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + abs(offset);

us = us * geo.DSO(1)/geo.DSD(1);
abstheta = abs(theta * geo.DSO(1)/geo.DSD(1));

w = ones(size(proj(:,:,1)));

for ii = 1:geo.nDetector
    t = us(ii);
    if(abs(t) <= abstheta)
        w(:,ii) = 0.5*(sin((pi/2)*atan(t/geo.DSO(1))/(atan(abstheta/geo.DSO(1)))) + 1);
    end
    if(t<-abstheta)
        w(:,ii) = 0;
    end
end

if(theta<0)
    w = fliplr(w);
end
proj_w=proj;% preallocation
for ii = 1:size(proj,3)
    proj_w(:,:,ii) = proj(:,:,ii).*w*2;
end

end

function [filter, parker,wang, gpuids]=parse_inputs(proj,geo,angles,argin)

opts =  {'filter','parker','wang', 'gpuids'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:FDK:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
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
            error('TIGRE:FDK:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
        end
        val=argin{jj};
    end
    
    switch opt
        case 'parker'
            if default
                if size(angles,1)==1 || (all(angles(2,:)==0) && all(angles(3,:)==0))
                    parker=max(angles(1,:))-min(angles(1,:))<(2*pi-max(diff(angles(1,:))));
                else
                    parker=false;
                end
            else
                parker=val;
            end
        case 'wang'
            if default
                wang=apply_wang_weights(geo);
            else
                wang=val;
            end
                    
  
        case 'filter'
            if default
                filter='ram-lak';
            else
                if  ~ischar( val)
                    error('TIGRE:FDK:InvalidInput','Invalid filter')
                end
                filter=val;
            end
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
       
        otherwise
            error('TIGRE:FDK:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in FDK()']);
    end
end
end


function bool = apply_wang_weights(geo)
    if (size(geo.offDetector,2) > 1) && length(unique(geo.offDetector(1,:)))>1
        warning('FDK Wang weights: varying offDetector detected, Wang weights not being applied');
        bool = false;
        return
    end
    
    if geo.offDetector(1) == 0
        bool = false;
        return
    end
    
    if (numel(geo.DSO) > 1) && (length(unique(geo.DSO))>1)
        warning('FDK Wang weights: varying DSO detected, Wang weights not being applied');
        bool = false;
        return
    end

    percent_offset = abs(geo.offDetector(1)/geo.sDetector(1)) * 100;    
    if percent_offset > 30
        warning('FDK Wang weights: Detector offset percent: %0.2f is greater than 30 which may result in image artifacts, consider rebinning 360 degree projections to 180 degrees', percent_offset)
    end
    
    bool = true;
end
