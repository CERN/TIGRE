function [image,residual,quality]=OS_SART_surface(data,geo,surface,niter,varargin)
% Explicit opt-in OS_SART; conventional calls retain upstream behavior.
[image,residual,quality]=allSurfaceRecon('OS_SART',data,geo,surface,niter,varargin{:});
end
