function [image,residual,quality]=hybrid_fLSQR_TV_surface(data,geo,surface,niter,varargin)
% Explicit opt-in hybrid_fLSQR_TV; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('hybrid_fLSQR_TV',data,geo,surface,niter,varargin{:});
end
