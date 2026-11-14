function [image,residual,quality]=hybrid_LSQR_surface(data,geo,surface,niter,varargin)
% Explicit opt-in hybrid_LSQR; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('hybrid_LSQR',data,geo,surface,niter,varargin{:});
end
