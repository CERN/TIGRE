function [image,residual,quality]=IRN_TV_CGLS_surface(data,geo,surface,niter,varargin)
% Explicit opt-in IRN_TV_CGLS; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('IRN_TV_CGLS',data,geo,surface,niter,varargin{:});
end
