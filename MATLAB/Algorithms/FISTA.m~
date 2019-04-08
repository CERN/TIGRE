function [res] = FISTA(proj,geo,angles,niter,hyper)
% FISTA is a quadratically converging algorithm that relies on the hyper
% parameter named 'hyper'. This parameter should approximate the largest 
% eigenvalue in the A matrix in the equation Ax-b and Atb. Empirical tests
% show that for, the headphantom object:
%           geo.nVoxel = [64,64,64]'    ,      hyper (approx=) 2.e8
%           geo.nVoxel = [512,512,512]' ,      hyper (approx=) 2.e4
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

lambda = 0.1;
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
    x_rec = im3DDenoise(res,'TV',20,1/lambdaforTV);  
    told = t;
    t = ( 1+sqrt(1+4*t*t) ) / 2;
    res= x_rec + (told-1)/t * (x_rec - x_recold);
    if (ii==1);
        expected_time=toc*niter;
        disp('FISTA');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
end

end
