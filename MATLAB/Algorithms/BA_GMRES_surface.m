function [image,residual,quality]=BA_GMRES_surface(data,geo,surface,niter,varargin)
% Explicit opt-in BA_GMRES; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('BA_GMRES',data,geo,surface,niter,varargin{:});
end
