function [image,residual,quality]=ISTA_surface(data,geo,surface,niter,varargin)
% Explicit opt-in ISTA; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('ISTA',data,geo,surface,niter,varargin{:});
end
