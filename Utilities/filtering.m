function [ proj ] = filtering(proj,geo,alpha)
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
%                     https://github.com/CERN/TIGRE/license.txt
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Kyungsang Kim, modified by Ander Biguri 
%--------------------------------------------------------------------------




filt_len = max(64,2^nextpow2(2*geo.nDetector(1)));
[ramp_kernel] = ramp_flat(filt_len);

d = 1; % cut off (0~1)
[filt] = Filter(geo.filter, ramp_kernel, filt_len, d);
filt = repmat(filt',[1 geo.nDetector(2)]);

for ii=1:length(alpha)
    
    fproj = (zeros(filt_len,geo.nDetector(2),'single'));
    
    fproj(filt_len/2-geo.nDetector(1)/2+1:filt_len/2+geo.nDetector(1)/2,:) = proj(:,:,ii);
    
    fproj = fft(fproj);   
    
    fproj = fproj.*filt;
    
    fproj = (real(ifft(fproj)));
    
  
    proj(:,:,ii) = fproj(end/2-geo.nDetector(1)/2+1:end/2+geo.nDetector(1)/2,:)/2/geo.dDetector(1)*(2*pi/  length(alpha)   )/2*(geo.DSD/geo.DSO);
    
    
end

proj=permute(proj,[2 1 3]);
end

function [h, nn] = ramp_flat(n)
nn = [-(n/2):(n/2-1)]';
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