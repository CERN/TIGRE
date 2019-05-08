function [res]=FDK(proj,geo,angles,varargin)
%TODO docs FDK
% 
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
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Kyungsang Kim, modified by Ander Biguri 
%--------------------------------------------------------------------------
[filter,parker]=parse_inputs(proj,geo,angles,varargin);
geo=checkGeo(geo,angles);
geo.filter=filter;

%Input is data,geosize,angles



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
    [uu,vv] = meshgrid(us,vs); %detector
    
    %Create weight according to each detector element
    w = (geo.DSD(ii))./sqrt((geo.DSD(ii))^2+uu.^2 + vv.^2);
    
    %Multiply the weights with projection data
    proj(:,:,ii) = proj(:,:,ii).*w';
end
%% filter
proj = filtering(proj,geo,angles,parker); % Not sure if offsets are good in here
%RMFIELD Remove fields from a structure array.
geo=rmfield(geo,'filter');
%% backproject
res=Atb((proj),geo,angles); % Weighting is inside


end

function [filter, parker]=parse_inputs(proj,geo,alpha,argin)
opts=     {'filter','parker'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('CBCT:FDK:InvalidInput','Invalid number of inputs')
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
    % if one option isnot default, then extranc value from input
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
        % % % % % % % Verbose
        case 'parker'
            if default
                parker=0;
            else
                parker=val;
            end
           
        case 'filter'
            if default
                filter='ram-lak';
            else
                if  ~ischar( val)
                    error('CBCT:FDK:InvalidInput','Invalid filter')
                end
                filter=val;
            end
       
        otherwise
            error('CBCT:FDK:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in FAK()']);
    end
end

end