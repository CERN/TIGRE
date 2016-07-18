function [res]=MLEM(proj,geo,angles,niter,varargin)
%MLEM solves the tomographic problem by using Maximum Likelihood Expection
% Maximitation algorithm. 
%
%   MLEM(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
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
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------

res=ones(geo.nVoxel.','single');
W=Atb(ones([geo.nDetector.' numel(angles)],'single'),geo,angles);

for ii=1:niter
   
    auxMLEM=proj./Ax(res,geo,angles);
    auxMLEM(isnan(auxMLEM)) = 0;
    auxMLEM(isinf(auxMLEM)) = 0;
    
    imgupdate = Atb(auxMLEM, geo,angles)./W;
    imgupdate(isnan(imgupdate)) = 0;
    imgupdate(isinf(imgupdate)) = 0;
    res = max(res.*imgupdate,0);

end

end
