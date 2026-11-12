function [image,residual,quality]=AB_GMRES_surface(data,geo,surface,niter,varargin)
% Explicit opt-in AB_GMRES; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('AB_GMRES',data,geo,surface,niter,varargin{:});
end
