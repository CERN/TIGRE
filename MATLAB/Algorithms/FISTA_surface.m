function [image,residual,quality]=FISTA_surface(data,geo,surface,niter,varargin)
% Explicit opt-in FISTA; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('FISTA',data,geo,surface,niter,varargin{:});
end
