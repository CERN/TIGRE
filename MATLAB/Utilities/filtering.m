function [ proj ] = filtering(proj,geo,angles,parker,varargin)
%FILTERING Summary of this function goes here
%   Detailed explanation goes here
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
[usegpufft,gpuids]=parse_inputs(varargin);

if parker
	proj = ParkerWeight(proj,geo,angles,parker);
	diff_angles = diff(angles(1,:)); angle_step = mean(abs(diff_angles)); % to be used later
end 

filt_len = max(64,2^nextpow2(2*geo.nDetector(1)));
[ramp_kernel] = ramp_flat(filt_len);

d = 1; % cut off (0~1)
[filt] = Filter(geo.filter, ramp_kernel, filt_len, d);
if usegpufft>0
    filt = filt';
else
    filt = repmat(filt',[1 geo.nDetector(2)]);
end

proj=permute(proj,[2 1 3]);

if usegpufft==2
    bundle_size = 32;  %len(gpuids)
    n_bundles = floor((size(angles,2)+bundle_size-1)/bundle_size);
    n_angles = size(angles, 2);
    for ii=1:n_bundles
        
        if ii ~= n_bundles
            bundle_size_actual = bundle_size;
        else
            bundle_size_actual = n_angles-(ii-1)*bundle_size;
        end
        idx_begin = (ii-1)*bundle_size+1;
        idx_end = (ii-1)*bundle_size+bundle_size_actual+1;
        proj_flt = ApplyPaddingAndFbpFiltration(proj, idx_begin, idx_end, filt, 1, gpuids.devices);
        bundle_range = idx_begin:(idx_end-1);
        if parker
            proj_flt = proj_flt/2/geo.dDetector(1)*(2*pi/  (pi/angle_step)  )/2*(geo.DSD(bundle_range)/geo.DSO(bundle_range));
        else
            proj_flt = proj_flt/2/geo.dDetector(1)*(2*pi/  size(angles,2)   )/2*(geo.DSD(bundle_range)/geo.DSO(bundle_range));
        end
        proj(:,:,bundle_range) = proj_flt;
    end
else 
    
    for ii=1:size(angles,2)

        fproj = (zeros(filt_len,geo.nDetector(2),'single'));

        fproj(round(filt_len/2-geo.nDetector(1)/2+1):round(filt_len/2+geo.nDetector(1)/2),:) = proj(:,:,ii);

        fproj = fft(fproj);   
        fproj = fproj.*filt;
        fproj = (real(ifft(fproj)));

        if parker
            proj(:,:,ii) = fproj(round(end/2-geo.nDetector(1)/2+1):round(end/2+geo.nDetector(1)/2),:)/2/geo.dDetector(1)*(2*pi/  (pi/angle_step)   )/2*(geo.DSD(ii)/geo.DSO(ii));
        else
            proj(:,:,ii) = fproj(round(end/2-geo.nDetector(1)/2+1):round(end/2+geo.nDetector(1)/2),:)/2/geo.dDetector(1)*(2*pi/  size(angles,2)   )/2*(geo.DSD(ii)/geo.DSO(ii));
        end 

    end
end

proj=permute(proj,[2 1 3]);
end

function [h, nn] = ramp_flat(n)
nn = (-(n/2):(n/2-1))';
h = zeros(size(nn),'single');
h(n/2+1) = 1 / 4;
odd = mod(nn,2) == 1;
h(odd) = -1 ./ (pi * nn(odd)).^2;
end


function [filt] = Filter(filter, kernel, order, d)

f_kernel = abs(fft(kernel))*2;
filt = f_kernel(1:order/2+1)';
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

switch lower(filter)
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'  
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        disp(filter);
        error('Invalid filter selected.');
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt , filt(end-1:-1:2)];    % Symmetry of the filter
return

end

function [usegpufft, gpuids]=parse_inputs(argin)

opts =  {'usegpufft','gpuids'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:filtering:InvalidInput','Invalid number of inputs')
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
            error('TIGRE:filtering:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
        end
        val=argin{jj};
    end
    
    switch opt
        case 'usegpufft'
            if default
                usegpufft=2;
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
            error('TIGRE:filtering:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in FDK()']);
    end
end

end
