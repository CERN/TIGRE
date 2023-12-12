function [sigma] = power_method(geo,angles)
tol=1e-6;
maxiter=100;
arr_old=randn(geo.nVoxel','single');
error = tol + 1;
i = 0;
while error >= tol

    % very verbose and inefficient for now
    omega = Ax(arr_old,geo,angles);
    alpha = im3Dnorm(omega,'L2');
    u = (1.0 / alpha) * omega;
    z = Atb(u,geo,angles,'matched');
    beta = im3Dnorm(z,'L2');
    arr = (1.0 / beta) * z;
    error = im3Dnorm(Ax(arr,geo,angles) - beta * u,'L2'); 
    sigma = beta;
    arr_old = arr;
    i = i+1;
    if i >= maxiter
        return
    end
end

end

