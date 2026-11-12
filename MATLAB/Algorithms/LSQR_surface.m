function [image,residual,quality]=LSQR_surface(data,geo,surface,niter,varargin)
% Explicit opt-in LSQR; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('LSQR',data,geo,surface,niter,varargin{:});
end
