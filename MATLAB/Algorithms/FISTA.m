function [res] = FISTA(proj,geo,angles,niter,hyper)
tic
lambda = 0.1;
res = zeros(geo.nVoxel','single');
x_rec = res;
L = hyper;
bm = 1/L;
t = 1;
for ii = 1:niter
    if (ii==1);tic;end
    % gradient descent step
    res = res + bm * 2 * Atb(proj - Ax(res,geo,angles, 'interpolated'), geo, angles, 'matched');
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
