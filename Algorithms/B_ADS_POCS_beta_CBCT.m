function [ fres ] = B_ADS_POCS_beta_CBCT(proj,geo,angles,iterPOCS,maxiter,epsilon,c)
% http://ieeexplore.ieee.org/xpl/abstractAuthors.jsp?arnumber=5874264
% ADS-POCS with Bregman Iteration

beta=1;
proj_k=proj;
for ii=1:maxiter
    [fres]= ADS_POCS_CBCT(proj_k,geo,angles,iterPOCS,epsilon);
    proj_k=proj_k+beta*c*(proj-Ax(fres,geo,angles,'Krylov'));
end