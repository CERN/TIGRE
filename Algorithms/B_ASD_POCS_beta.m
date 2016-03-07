function [ fres ] = B_ASD_POCS_beta(proj,geo,angles,iterPOCS,maxiter,epsilon,c)
% http://ieeexplore.ieee.org/xpl/abstractAuthors.jsp?arnumber=5874264
% ADS-POCS with Bregman Iteration

beta=1;
proj_k=proj;
for ii=1:maxiter
    [fres]= ADS_POCS_CBCT(proj_k,geo,angles,iterPOCS,epsilon);
    proj_k=proj_k+beta*c*(proj-Ax(fres,geo,angles,'Krylov'));
end