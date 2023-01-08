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
[filter,parker,dowang,usegpufft,gpuids]=parse_inputs(angles,varargin);

geo=checkGeo(geo,angles);
geo.filter=filter;

if dowang
    % Zero-padding to avoid FFT-induced aliasing %TODO: should't this be
    % for all cases, not just wang?
    % Preweighting using Wang function
    proj=proj.*redundancy_weighting(geo);
    [proj, geo] = zeropadding(proj, geo);

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
proj = filtering(proj,geo,angles,parker,usegpufft); % Not sure if offsets are good in here

% RMFIELD Remove fields from a structure array.
geo=rmfield(geo,'filter');
%% backproject
%%%
% [proj, w] = preweighting(proj,geo);
% imshow(w,[])
%%%
res=Atb((proj),geo,angles, 'gpuids', gpuids); % Weighting is inside


end

function [zproj, zgeo] = zeropadding(proj, geo)
% ZEROPADDING as preprocessing for preweighting
zgeo = geo;

offset = geo.offDetector(1);
offset = offset + (geo.DSD(1) / geo.DSO(1)) * geo.COR(1);

padwidth = fix(2*offset./geo.dDetector(1));
zgeo.offDetector(1,:) = geo.offDetector(1,:) - padwidth/2 * geo.dDetector(1);
zgeo.nDetector(1) = abs(padwidth) + geo.nDetector(1);
zgeo.sDetector(1) = zgeo.nDetector(1) * zgeo.dDetector(1);

% Pad on the left size when offset >0
if(offset>0)
    for ii = 1:size(proj,3)
        zproj(:,:,ii) = [zeros(size(proj,1), padwidth), proj(:,:,ii)];
    end
else
    for ii = 1:size(proj,3)
        zproj(:,:,ii) = [proj(:,:,ii), zeros(size(proj,1), abs(padwidth))];
    end
end

end


function [filter, parker, wang, usegpufft, gpuids]=parse_inputs(angles,argin)

opts =  {'filter','parker','wang','usegpufft','gpuids'};
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
            error('CBCT:FDK:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
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
                wang=true;
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
        case 'usegpufft'
            if default
                usegpufft=true;
            else
                usegpufft=val;
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
